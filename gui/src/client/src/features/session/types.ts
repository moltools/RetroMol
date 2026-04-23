import { z } from "zod";

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

export function initSession(): Session {
  const newSession = SessionSchema.parse({});
  return newSession;
};