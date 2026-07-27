import { postJson } from "../http";
import { Session, initSession, SessionSchema, CreateSessionRespSchema, GetSessionRespSchema, SseTicketRespSchema } from "./types";
import { getCookie } from "./utils";
import { z } from "zod";

const OkRespSchema = z.object({ ok: z.boolean() }).partial();

export async function createSession(): Promise<string> {
  const newSession = initSession();
  const data = await postJson("/api/createSession", { session: newSession }, CreateSessionRespSchema);
  return data.sessionId;
};

export async function getSession(sessionIdArg?: string): Promise<Session> {
  const sessionId = sessionIdArg ?? getCookie("sessionId");
  if (!sessionId) throw new Error("No sessionId provided or found in cookies");
  const data = await postJson("/api/getSession", { sessionId }, GetSessionRespSchema);
  // Dont have to sanitize session: postJson already validates with GetSessionRespSchema
  return data.session;
};

export async function refreshSession(sessionId: string): Promise<Session> {
  return getSession(sessionId);
};

export async function saveSession(session: Session): Promise<void> {
  // Runtime validate before sending (especially useful because session is user-mutated in UI)
  const sanitized = SessionSchema.parse(session);
  await postJson("/api/saveSession", { session: sanitized }, z.unknown());
};

export async function deleteSession(): Promise<void> {
  const sessionId = getCookie("sessionId");
  if (!sessionId) throw new Error("No sessionId cookie found");
  await postJson("/api/deleteSession", { sessionId }, z.unknown());
};

export async function deleteSessionItem(sessionId: string, itemId: string): Promise<void> {
  await postJson("/api/deleteSessionItem", { sessionId, itemId }, OkRespSchema);
};

// EventSource can only issue a plain GET (no custom headers/body), so the SSE
// endpoint can't be authorized with a POSTed sessionId like everything else
// here. Instead we mint a short-lived, single-use ticket via a normal POST,
// and only that opaque ticket ever ends up in the SSE URL.
export async function getSessionEventsTicket(sessionId: string): Promise<string> {
  const data = await postJson("/api/sessionEventsTicket", { sessionId }, SseTicketRespSchema);
  return data.ticket;
};