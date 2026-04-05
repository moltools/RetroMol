import { defaultMotifColorMap } from "./utils";
import { z } from "zod";

export const BaseItemSchema = z.object({
  id: z.string(),
  name: z.string(), // display name
  status: z.enum(["queued", "processing", "done", "error"]).default("queued"),
  errorMessage: z.string().nullable().optional(),
  updatedAt: z.number().nonnegative().default(() => Date.now()),
})

export const PrimarySequenceMotifSchema = z.object({
  id: z.string(),
  name: z.string().nullable().optional(),
  displayName: z.string().nullable().optional(),

  // Structural information
  tags: z.array(z.number().int().nonnegative()).default([]),
  smiles: z.string().nullable().optional(),
  morganfingerprint2048r2: z.string().length(512).nullable().optional(),
}).refine(
  ({ smiles, morganfingerprint2048r2 }) =>
    (smiles == null && morganfingerprint2048r2 == null) ||
    (smiles != null && morganfingerprint2048r2 != null),
    { message: "Both 'smiles' and 'morganfingerprint2048r2' must be provided together or both be null",}
)

export type PrimarySequenceMotif = z.output<typeof PrimarySequenceMotifSchema>;

export const PrimarySequenceSchema = z.object({
  id: z.string(),
  name: z.string().nullable().optional(),
  parentSmilesTagged: z.string().nullable().optional(),
  sequence: z.array(PrimarySequenceMotifSchema).min(1),
})

export type PrimarySequence = z.output<typeof PrimarySequenceSchema>;

export const MsaSequenceSchema = PrimarySequenceSchema.extend({
  itemId: z.string(),
  primarySequenceId: z.string(),
  hidden: z.boolean().default(false),
});

export type MsaSequence = z.output<typeof MsaSequenceSchema>;

export const MsaStateSchema = z.object({
  aligned: z.boolean().default(false),
  centerId: z.string().nullable().optional(),
  sequences: z.array(MsaSequenceSchema).default([]),
});

export type MsaState = z.output<typeof MsaStateSchema>;

export const BaseFingerprintSchema = z.object({
  id: z.string(),
  retrofingerprint512: z.string().length(128),
  score: z.number().min(0).max(1),
})

export const CompoundRetrofingerprintSchema = BaseFingerprintSchema.extend({});
export const GeneClusterRetrofingerprintSchema = BaseFingerprintSchema.extend({});

export const CompoundItemSchema = BaseItemSchema.extend({
  kind: z.literal("compound"),
  smiles: z.string(),
  taggedSmiles: z.string().nullable().optional(),
  retrofingerprints: z.array(CompoundRetrofingerprintSchema).default([]),
  primarySequences: z.array(PrimarySequenceSchema).default([]),
})

export const GeneClusterSchema = BaseItemSchema.extend({
  kind: z.literal("gene_cluster"),
  fileContent: z.string(),
  retrofingerprints: z.array(GeneClusterRetrofingerprintSchema).default([]),
  primarySequences: z.array(PrimarySequenceSchema).default([]),
})

export const SessionItemSchema = z.discriminatedUnion("kind", [
  CompoundItemSchema,
  GeneClusterSchema,
])

export type CompoundItem = z.output<typeof CompoundItemSchema>;
export type GeneClusterItem = z.output<typeof GeneClusterSchema>;
export type SessionItem = z.output<typeof SessionItemSchema>;

export const AlignmentTypeSchema = z.enum(["global", "local"]);
export const EmbeddingVisualizationTypeSchema = z.enum(["pca", "umap"]);
export const QuerySearchSpaceSchema = z.enum(["only_compounds", "only_gene_clusters", "both"]);
export const AnnotationFilterSchema = z.object({scheme: z.string(), key: z.string(),  value: z.string()});

export type AlignmentType = z.output<typeof AlignmentTypeSchema>;
export type EmbeddingVisualizationType = z.output<typeof EmbeddingVisualizationTypeSchema>;
export type QuerySearchSpace = z.output<typeof QuerySearchSpaceSchema>;
export type AnnotationFilter = z.output<typeof AnnotationFilterSchema>;

export const MsaSettingsSchema = z.object({
  alignmentType: AlignmentTypeSchema.default("global"),
});

export type MsaSettings = z.output<typeof MsaSettingsSchema>;

export const QuerySettingsSchema = z.object({
  similarityThreshold: z.number().min(0).max(1).default(0.7),
  searchSpace: QuerySearchSpaceSchema.default("both"),
  annotationFilters: z.array(AnnotationFilterSchema).default([]),
});

export type QuerySettings = z.output<typeof QuerySettingsSchema>;

export const SessionSettingsSchema = z.object({
  motifColorPalette: z.record(z.string()).default(() => defaultMotifColorMap()),
  embeddingVisualizationType: EmbeddingVisualizationTypeSchema.default("pca"),
  msaSettings: MsaSettingsSchema.default(() => ({})),
  querySettings: QuerySettingsSchema.default(() => ({})),
});

export type SessionSettings = z.output<typeof SessionSettingsSchema>;

export const SessionSchema = z.object({
  sessionId: z.string().default(() => crypto.randomUUID()),
  created: z.number().nonnegative().default(() => Date.now()),
  items: z.array(SessionItemSchema).default([]),
  settings: SessionSettingsSchema.default(() => ({})),
  msaState: MsaStateSchema.default(() => ({})),
})

export type Session = z.output<typeof SessionSchema>;

// Simple response wrappers
export const CreateSessionRespSchema = z.object({ sessionId: z.string() });
export const GetSessionRespSchema = z.object({ session: SessionSchema });

export function initSession(): Session {
  const newSession = SessionSchema.parse({});
  return newSession;
}