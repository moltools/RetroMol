import { z } from "zod";
import { PrimarySequenceItemSchema } from "../reconstruction/types";

export const BaseItemSchema = z.object({
  id: z.string(),
  name: z.string(), // display name
  score: z.number().min(0).max(1).nullable().optional(),
  status: z.enum(["queued", "processing", "done", "error"]).default("queued"),
  errorMessage: z.string().nullable().optional(),
  updatedAt: z.number().nonnegative().default(() => Date.now()),
});

export const CompoundItemSchema = BaseItemSchema.extend({
  kind: z.literal("compound"),
  smiles: z.string(),
  matchStereochemistry: z.boolean(),
  // User overrides for the primary sequence(s) RetroMol parsed out of this
  // compound, keyed by reconstruction index (stringified, for JSON-safe object
  // keys). A missing key means "use the algorithm's own parse" -- reverting an
  // edit just deletes that key rather than storing a copy of the original, so
  // the original can never go stale.
  editedPrimarySequences: z.record(z.string(), z.array(PrimarySequenceItemSchema)).nullable().optional(),
});

export const ClusterItemSchema = BaseItemSchema.extend({
  kind: z.literal("cluster"),
  fileContent: z.string(),
});

export const SessionItemSchema = z.discriminatedUnion("kind", [CompoundItemSchema, ClusterItemSchema]);

export type CompoundItem = z.output<typeof CompoundItemSchema>;
export type ClusterItem = z.output<typeof ClusterItemSchema>;
export type SessionItem = z.output<typeof SessionItemSchema>;

export const SessionSchema = z.object({
  sessionId: z.string().default(() => crypto.randomUUID()),
  created: z.number().nonnegative().default(() => Date.now()),
  items: z.array(SessionItemSchema).default([]),
});

export type Session = z.output<typeof SessionSchema>;

// Simple response wrappers
export const CreateSessionRespSchema = z.object({ sessionId: z.string() });
export const GetSessionRespSchema = z.object({ session: SessionSchema });
export const SseTicketRespSchema = z.object({
  ticket: z.string(),
  expiresInSeconds: z.number().optional(),
});

export function initSession(): Session {
  const newSession = SessionSchema.parse({});
  return newSession;
};