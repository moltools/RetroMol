import { useEffect, useState } from "react";

type Listener = () => void;

const listeners = new Set<Listener>();
let intervalId: number | null = null;

function ensureInterval(periodMs: number) {
  if (intervalId !== null) return;
  intervalId = window.setInterval(() => {
    listeners.forEach((listener) => listener());
  }, periodMs);
}

function maybeStopInterval() {
  if (listeners.size === 0 && intervalId !== null) {
    window.clearInterval(intervalId);
    intervalId = null;
  }
}

/**
 * Subscribes to a single shared interval timer, rather than each caller
 * starting its own setInterval. Meant for "updated 3s ago"-style UI that
 * needs to re-render periodically without any real state change.
 *
 * Assumes all concurrent callers want (roughly) the same period — the period
 * used is whichever caller happens to start the shared timer first.
 */
export function useTick(periodMs: number = 5000): void {
  const [, setTick] = useState(0);

  useEffect(() => {
    ensureInterval(periodMs);
    const listener = () => setTick((t) => t + 1);
    listeners.add(listener);
    return () => {
      listeners.delete(listener);
      maybeStopInterval();
    };
  }, [periodMs]);
}
