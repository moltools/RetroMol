import React from "react";
import { act, renderHook } from "@testing-library/react";
import { deleteSessionItem, refreshSession, saveSession } from "./api";
import { postJson } from "../http";
import type { Session } from "./types";
import { useSessionPolling } from "./useSessionPolling";

jest.mock("./api", () => ({ ...jest.requireActual("./api"), refreshSession: jest.fn() }));
jest.mock("../http");
const refresh = refreshSession as jest.Mock;
const post = postJson as jest.Mock;

function makeSession(status: "queued" | "processing" | "done" = "queued", sessionId = "s1"): Session {
  return {
    sessionId, created: 1,
    items: [{ id: "i1", name: "Compound", kind: "compound", smiles: "C", matchStereochemistry: false, status, updatedAt: 1 }],
  };
}

function setup(initial = makeSession()) {
  const onError = jest.fn();
  const hook = renderHook(() => {
    const [session, setSession] = React.useState<Session | null>(initial);
    useSessionPolling(session, setSession, onError);
    return { session, setSession };
  });
  return { ...hook, onError };
}

async function advance(ms: number) {
  await act(async () => { jest.advanceTimersByTime(ms); });
}

async function setHidden(hidden: boolean) {
  await act(async () => {
    Object.defineProperty(document, "hidden", { configurable: true, value: hidden });
    document.dispatchEvent(new Event("visibilitychange"));
  });
}

describe("workspace session polling", () => {
  beforeEach(() => {
    jest.useFakeTimers();
    jest.resetAllMocks();
    Object.defineProperty(document, "hidden", { configurable: true, value: false });
    refresh.mockResolvedValue(makeSession());
  });

  afterEach(() => {
    jest.useRealTimers();
    delete (document as any).hidden;
  });

  it("refreshes job progress, then slows down when all jobs finish", async () => {
    refresh.mockResolvedValueOnce(makeSession("processing")).mockResolvedValue(makeSession("done"));
    const { result } = setup();
    await advance(2_000);
    expect(result.current.session?.items[0].status).toBe("processing");
    await advance(2_000);
    expect(result.current.session?.items[0].status).toBe("done");
    await advance(9_999);
    expect(refresh).toHaveBeenCalledTimes(2);
    await advance(1);
    expect(refresh).toHaveBeenCalledTimes(3);
  });

  it("never overlaps requests when a refresh is slow", async () => {
    let resolve!: (session: Session) => void;
    refresh.mockReturnValueOnce(new Promise<Session>((done) => { resolve = done; }));
    setup();
    await advance(2_000);
    await advance(8_000);
    expect(refresh).toHaveBeenCalledTimes(1);
    await act(async () => { resolve(makeSession()); });
    await advance(2_000);
    expect(refresh).toHaveBeenCalledTimes(2);
  });

  it("pauses hidden tabs and refreshes immediately when visible again", async () => {
    setup(makeSession("done"));
    await setHidden(true);
    await advance(60_000);
    expect(refresh).not.toHaveBeenCalled();
    await setHidden(false);
    expect(refresh).toHaveBeenCalledTimes(1);
  });

  it("aborts on unmount and does not schedule another request", async () => {
    let resolve!: (session: Session) => void;
    refresh.mockReturnValueOnce(new Promise<Session>((done) => { resolve = done; }));
    const { unmount, onError } = setup();
    await advance(2_000);
    const signal = refresh.mock.calls[0][1] as AbortSignal;
    unmount();
    expect(signal.aborted).toBe(true);
    await act(async () => { resolve(makeSession()); });
    await advance(60_000);
    expect(refresh).toHaveBeenCalledTimes(1);
    expect(onError).not.toHaveBeenCalled();
  });

  it("does not overwrite local changes made while a refresh is in flight", async () => {
    let resolve!: (session: Session) => void;
    refresh.mockReturnValueOnce(new Promise<Session>((done) => { resolve = done; }));
    const { result } = setup();
    await advance(2_000);
    const local = makeSession();
    local.items[0].name = "Edited locally";
    act(() => result.current.setSession(local));
    await act(async () => { resolve(makeSession()); });
    expect(result.current.session).toBe(local);
  });

  it("discards late responses after switching sessions", async () => {
    let resolve!: (session: Session) => void;
    refresh.mockReturnValueOnce(new Promise<Session>((done) => { resolve = done; }));
    const { result } = setup();
    await advance(2_000);
    const signal = refresh.mock.calls[0][1] as AbortSignal;
    const next = makeSession("queued", "s2");
    act(() => result.current.setSession(next));
    expect(signal.aborted).toBe(true);
    await act(async () => { resolve(makeSession()); });
    expect(result.current.session).toBe(next);
    await advance(2_000);
    expect(refresh.mock.calls[1][0]).toBe("s2");
  });

  it("waits for an optimistic save to finish before polling", async () => {
    let finishSave!: () => void;
    post.mockReturnValueOnce(new Promise<void>((resolve) => { finishSave = resolve; }));
    setup();
    const saving = saveSession(makeSession());
    await advance(6_000);
    expect(refresh).not.toHaveBeenCalled();
    await act(async () => { finishSave(); await saving; });
    await advance(2_000);
    expect(refresh).toHaveBeenCalledTimes(1);
  });

  it("discards a snapshot if a server write completed while it was in flight", async () => {
    let resolve!: (session: Session) => void;
    refresh.mockReturnValueOnce(new Promise<Session>((done) => { resolve = done; }));
    const { result } = setup();
    const before = result.current.session;
    await advance(2_000);
    // No React state change here: the write-version guard must catch this race.
    await deleteSessionItem("s1", "i1");
    await act(async () => { resolve(makeSession("done")); });
    expect(result.current.session).toBe(before);
  });

  it("keeps polling paused until all concurrent writes settle, including failures", async () => {
    let finishSave!: () => void;
    let failDelete!: (error: Error) => void;
    post.mockReturnValueOnce(new Promise<void>((resolve) => { finishSave = resolve; }))
      .mockReturnValueOnce(new Promise<void>((_resolve, reject) => { failDelete = reject; }));
    setup();
    const saving = saveSession(makeSession());
    const deleting = deleteSessionItem("s1", "i1").catch(() => undefined);
    await act(async () => { finishSave(); await saving; });
    await advance(2_000);
    expect(refresh).not.toHaveBeenCalled();
    await act(async () => { failDelete(new Error("Delete failed")); await deleting; });
    await advance(2_000);
    expect(refresh).toHaveBeenCalledTimes(1);
  });

  it("backs off on failures, reports once, and resumes normal polling after recovery", async () => {
    refresh.mockRejectedValue(new Error("Server unavailable"));
    const { onError } = setup();
    await advance(2_000);
    await advance(3_999);
    expect(refresh).toHaveBeenCalledTimes(1);
    await advance(1);
    expect(refresh).toHaveBeenCalledTimes(2);
    await advance(8_000);
    expect(refresh).toHaveBeenCalledTimes(3);
    await advance(16_000);
    expect(refresh).toHaveBeenCalledTimes(4);
    await advance(29_999);
    expect(refresh).toHaveBeenCalledTimes(4);
    expect(onError).toHaveBeenCalledTimes(1);
    refresh.mockResolvedValue(makeSession());
    await advance(1);
    await advance(2_000);
    expect(refresh).toHaveBeenCalledTimes(6);
  });

  it("times out a stalled request and retries with backoff", async () => {
    refresh.mockImplementationOnce((_id, signal: AbortSignal) => new Promise((_resolve, reject) => {
      signal.addEventListener("abort", () => reject(new Error("Aborted")));
    }));
    const { onError } = setup();
    await advance(2_000);
    await advance(10_000);
    expect((refresh.mock.calls[0][1] as AbortSignal).aborted).toBe(true);
    expect(onError).toHaveBeenCalledWith(expect.stringContaining("Request timed out"));
    await advance(4_000);
    expect(refresh).toHaveBeenCalledTimes(2);
  });
});
