import { postJson } from "../http";
import { Session, initSession, SessionSchema, CreateSessionRespSchema, GetSessionRespSchema } from "./types";
import { getCookie } from "./utils";
import { z } from "zod";

export async function createSession(): Promise<string> {
  const newSession = initSession();
  const data = await postJson("/api/createSession", { session: newSession }, CreateSessionRespSchema);
  return data.sessionId;
}

export async function getSession(sessionIdArg?: string): Promise<Session> {
  const sessionId = sessionIdArg ?? getCookie("sessionId");
  if (!sessionId) throw new Error("No sessionId provided or found in cookies");
  const data = await postJson("/api/getSession", { sessionId }, GetSessionRespSchema);
  return data.session;
}

export async function refreshSession(sessionId: string): Promise<Session> {
  return getSession(sessionId);
}

export async function saveSession(session: Session): Promise<void> {
  // Runtime validate before sending (especially useful because session is user-mutated in UI)
  SessionSchema.parse(session);
  await postJson("/api/saveSession", { session }, z.unknown());
}

export async function deleteSession(): Promise<void> {
  const sessionId = getCookie("sessionId");
  if (!sessionId) throw new Error("No sessionId cookie found");
  await postJson("/api/deleteSession", { sessionId }, z.unknown());
}