import { postJson } from "../http";
import { Session, initSession, SessionSchema, CreateSessionRespSchema, GetSessionRespSchema } from "./types";
import { getCookie } from "./utils";
import { z } from "zod";

const OkRespSchema = z.object({ ok: z.boolean() }).partial();

let pendingSessionWrites = 0;
let sessionWriteVersion = 0;

// A poll may only apply its snapshot if no write was pending when it started
// and no write has begun since. Track writes independently of React renders,
// including optimistic edits whose save has not reached the server yet.
export function getSessionWriteVersion(): number | null {
  return pendingSessionWrites ? null : sessionWriteVersion;
}

export async function withSessionWrite<T>(operation: () => Promise<T>): Promise<T> {
  pendingSessionWrites += 1;
  sessionWriteVersion += 1;
  try {
    return await operation();
  } finally {
    pendingSessionWrites -= 1;
  }
}

export async function createSession(): Promise<string> {
  const newSession = initSession();
  const data = await postJson("/api/createSession", { session: newSession }, CreateSessionRespSchema);
  return data.sessionId;
};

export async function getSession(sessionIdArg?: string, signal?: AbortSignal): Promise<Session> {
  const sessionId = sessionIdArg ?? getCookie("sessionId");
  if (!sessionId) throw new Error("No sessionId provided or found in cookies");
  const data = await postJson("/api/getSession", { sessionId }, GetSessionRespSchema, signal);
  // Dont have to sanitize session: postJson already validates with GetSessionRespSchema
  return data.session;
};

export async function refreshSession(sessionId: string, signal?: AbortSignal): Promise<Session> {
  return getSession(sessionId, signal);
};

export async function saveSession(session: Session): Promise<void> {
  // Runtime validate before sending (especially useful because session is user-mutated in UI)
  const sanitized = SessionSchema.parse(session);
  await withSessionWrite(() => postJson("/api/saveSession", { session: sanitized }, z.unknown()));
};

export async function deleteSession(): Promise<void> {
  const sessionId = getCookie("sessionId");
  if (!sessionId) throw new Error("No sessionId cookie found");
  await withSessionWrite(() => postJson("/api/deleteSession", { sessionId }, z.unknown()));
};

export async function deleteSessionItem(sessionId: string, itemId: string): Promise<void> {
  await withSessionWrite(() => postJson("/api/deleteSessionItem", { sessionId, itemId }, OkRespSchema));
};

export async function renameSessionItem(session: Session, itemId: string, name: string): Promise<Session> {
  const nextSession: Session = {
    ...session,
    items: session.items.map((it) => (it.id === itemId ? { ...it, name } : it)),
  };
  await saveSession(nextSession);
  return nextSession;
};
