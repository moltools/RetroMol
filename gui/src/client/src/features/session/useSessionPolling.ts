import { useEffect, useRef } from "react";
import type { Dispatch, SetStateAction } from "react";
import { getSessionWriteVersion, refreshSession } from "./api";
import type { Session } from "./types";

// Ordinary requests release the Gunicorn thread between refreshes. A persistent
// SSE request occupied one of only eight API threads for every open workspace.
const BUSY_INTERVAL_MS = 2_000;
const IDLE_INTERVAL_MS = 10_000;
const MAX_BACKOFF_MS = 30_000;
const REQUEST_TIMEOUT_MS = 10_000;

export function useSessionPolling(
  session: Session | null,
  setSession: Dispatch<SetStateAction<Session | null>>,
  onError: (message: string) => void,
) {
  const latestSession = useRef(session);
  latestSession.current = session;
  const sessionId = session?.sessionId;
  const busy = session?.items.some((item) => item.status === "queued" || item.status === "processing") ?? false;

  useEffect(() => {
    if (!sessionId) return;

    let alive = true;
    let timer: number | undefined;
    let controller: AbortController | null = null;
    let failures = 0;
    const interval = busy ? BUSY_INTERVAL_MS : IDLE_INTERVAL_MS;

    const schedule = (delay: number) => {
      window.clearTimeout(timer);
      if (alive && !document.hidden) timer = window.setTimeout(poll, delay);
    };

    const poll = async () => {
      // Schedule only after the previous request settles: a slow server must
      // never build up overlapping polls from one workspace.
      if (!alive || document.hidden || controller) return;
      window.clearTimeout(timer);
      const writeVersion = getSessionWriteVersion();
      if (writeVersion === null) {
        schedule(interval);
        return;
      }
      const snapshot = latestSession.current;
      const request = new AbortController();
      controller = request;
      let timedOut = false;
      const timeout = window.setTimeout(() => {
        timedOut = true;
        request.abort();
      }, REQUEST_TIMEOUT_MS);

      try {
        const fresh = await refreshSession(sessionId, request.signal);
        if (!alive || request.signal.aborted) return;
        failures = 0;
        // Do not overwrite a local edit or submission made during this request.
        setSession((current) =>
          alive && current === snapshot && getSessionWriteVersion() === writeVersion ? fresh : current
        );
      } catch (err) {
        if (!alive || document.hidden || (request.signal.aborted && !timedOut)) return;
        failures += 1;
        if (failures === 1) {
          const message = request.signal.aborted
            ? "Request timed out"
            : err instanceof Error ? err.message : String(err);
          onError(`Failed to refresh session: ${message}. Retrying automatically.`);
        }
      } finally {
        window.clearTimeout(timeout);
        controller = null;
        schedule(failures ? Math.min(interval * 2 ** Math.min(failures, 5), MAX_BACKOFF_MS) : interval);
      }
    };

    const onVisibilityChange = () => {
      window.clearTimeout(timer);
      if (document.hidden) controller?.abort();
      else void poll();
    };

    schedule(interval);
    document.addEventListener("visibilitychange", onVisibilityChange);
    return () => {
      alive = false;
      window.clearTimeout(timer);
      controller?.abort();
      document.removeEventListener("visibilitychange", onVisibilityChange);
    };
  }, [sessionId, busy, setSession, onError]);
}
