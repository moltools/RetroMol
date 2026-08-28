import { z } from "zod";
import { PrimarySequenceItemSchema } from "../reconstruction/types";
import {
  DiscoveryEntryTypeSchema,
  DiscoveryMsaRowSchema,
  DiscoveryResultSchema,
  DiscoveryScoreModeSchema,
} from "../discovery/types";

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
  // Minimum PARAS prediction probability required to call a substrate for an NRPS
  // A-domain module (see retromol_antismash.inference.model_paras.ParasModel).
  // Doesn't affect PKS modules, which are classified directly from domain anatomy.
  parasThreshold: z.number().min(0).max(1).default(0.1),
  // User overrides for the primary sequence RetroMol parsed out of this gene
  // cluster, keyed by region index (stringified, for JSON-safe object keys) --
  // same convention as CompoundItemSchema.editedPrimarySequences above (a
  // missing key means "use the algorithm's own parse"; reverting just deletes
  // the key rather than storing a copy of the original).
  editedPrimarySequences: z.record(z.string(), z.array(PrimarySequenceItemSchema)).nullable().optional(),
});

export const DiscoveryQueryFlagsSchema = z.object({
  computeMsa: z.boolean().default(false),
  computeCompare: z.boolean().default(false),
  // Not precomputed (see run_discovery_query_job) -- gates whether the results
  // dialog's Enrichment view is enabled; the analysis itself is fetched live from
  // /api/enrichmentAnalysis based on whichever nearest neighbors are checked.
  computeEnrichment: z.boolean().default(true),
});
export type DiscoveryQueryFlags = z.output<typeof DiscoveryQueryFlagsSchema>;

export const DiscoveryQuerySettingsSchema = z.object({
  primarySequence: z.array(z.string()),
  entryType: DiscoveryEntryTypeSchema,
  scoreMode: DiscoveryScoreModeSchema,
  n: z.number().int(),
  topX: z.number().int(),
  includeUserUploads: z.boolean(),
  onlyUserUploads: z.boolean(),
  // The originating compound's own SMILES, captured at submit time -- see
  // WorkspaceDiscovery's queryOriginSmiles. Null when the query wasn't seeded from
  // a compound with a known structure (e.g. built from scratch, or from a cluster).
  queryOriginSmiles: z.string().nullable(),
});
export type DiscoveryQuerySettings = z.output<typeof DiscoveryQuerySettingsSchema>;

// Present once status === "done". Bundles the base ranking plus whichever optional
// metrics were flagged at submission (see DiscoveryQueryFlagsSchema) -- a metric
// missing here means it wasn't requested, not that it failed silently.
export const DiscoveryQueryPayloadSchema = z.object({
  querySequence: z.array(z.string()),
  selfScore: z.number(),
  candidatesConsidered: z.number().int().nonnegative(),
  candidatesSkipped: z.number().int().nonnegative(),
  results: z.array(DiscoveryResultSchema),
  msaRows: z.array(DiscoveryMsaRowSchema).optional(),
  compareValues: z.array(z.object({ id: z.string(), tanimotoSimilarity: z.number() })).optional(),
  compareRadius: z.number().optional(),
  compareNBits: z.number().optional(),
});
export type DiscoveryQueryPayload = z.output<typeof DiscoveryQueryPayloadSchema>;

export const DiscoveryQueryItemSchema = BaseItemSchema.extend({
  kind: z.literal("discoveryQuery"),
  settings: DiscoveryQuerySettingsSchema,
  flags: DiscoveryQueryFlagsSchema,
  payload: DiscoveryQueryPayloadSchema.nullable().optional(),
});
export type DiscoveryQueryItem = z.output<typeof DiscoveryQueryItemSchema>;

export const SessionItemSchema = z.discriminatedUnion("kind", [
  CompoundItemSchema,
  ClusterItemSchema,
  DiscoveryQueryItemSchema,
]);

export type CompoundItem = z.output<typeof CompoundItemSchema>;
export type ClusterItem = z.output<typeof ClusterItemSchema>;
export type SessionItem = z.output<typeof SessionItemSchema>;

export const SubmitDiscoveryQueryRespSchema = z.object({
  ok: z.boolean().optional(),
  item: DiscoveryQueryItemSchema,
});

// payload is stripped from every item by getSession (see routes/session.py) --
// this is the discoveryQuery equivalent of reconstructCompound/getClusterReadout,
// fetched on demand once a query's status is "done" (see DialogViewDiscoveryQuery).
export const GetDiscoveryQueryResultRespSchema = z.object({
  ok: z.boolean().optional(),
  status: BaseItemSchema.shape.status,
  payload: DiscoveryQueryPayloadSchema.nullable(),
});

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