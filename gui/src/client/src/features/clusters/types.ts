import { z } from "zod";

// Predicted substrate for an NRPS (A-domain) module, as returned by the PARAS model.
export const ClusterSubstrateSchema = z.object({
  name: z.string().nullable(),
  smiles: z.string().nullable(),
  score: z.number().nullable(),
});

export const ClusterNRPSAnatomySchema = z.object({
  has_C: z.boolean(),
  has_T: z.boolean(),
  has_E: z.boolean(),
  has_MT: z.boolean(),
  has_Ox: z.boolean(),
  has_R: z.boolean(),
  has_TE: z.boolean(),
});

export const ClusterPKSAnatomySchema = z.object({
  AT_loading_mode: z.enum(["cis", "trans", "unknown"]),
  has_active_KR: z.boolean(),
  has_active_DH: z.boolean(),
  has_active_ER: z.boolean(),
});

const ClusterModuleBaseSchema = z.object({
  module_index_in_gene: z.number(),
  start: z.number(),
  end: z.number(),
  gene_id: z.string(),
  gene_strand: z.enum(["+", "-"]),
  present_domains: z.array(z.string()),
});

export const ClusterModuleSchema = z.discriminatedUnion("type", [
  ClusterModuleBaseSchema.extend({
    type: z.literal("NRPS"),
    anatomy: ClusterNRPSAnatomySchema,
    predicted_substrate: ClusterSubstrateSchema.nullable().optional(),
  }),
  ClusterModuleBaseSchema.extend({
    type: z.literal("PKS"),
    anatomy: ClusterPKSAnatomySchema,
  }),
]);

// One antiSMASH region's linear module readout.
export const ClusterReadoutSchema = z.object({
  id: z.string(),
  file_name: z.string(),
  start: z.number(),
  end: z.number(),
  qualifiers: z.record(z.string(), z.any()).optional(),
  modules: z.array(ClusterModuleSchema),
  modifiers: z.array(z.string()),
});

export const ClusterPayloadSchema = z.object({
  readouts: z.array(ClusterReadoutSchema).default([]),
});

export const GetClusterReadoutRespSchema = z.object({
  ok: z.boolean().optional(),
  status: z.string().optional(),
  data: ClusterPayloadSchema,
});

// One antiSMASH region's readout, mapped onto the same matching-rule vocabulary a
// compound's primary sequence is drawn from (see
// retromol_antismash.modules.bgc_primary_sequence on the backend). Unlike a
// compound's reconstruction, this is just names -- no per-block tags, and no edited
// overrides yet.
export const ClusterPrimarySequenceSchema = z.object({
  id: z.string(),
  primary_sequence: z.array(z.string()),
});

export const ReconstructGeneClusterRespSchema = z.object({
  ok: z.boolean().optional(),
  status: z.string().optional(),
  data: z.array(ClusterPrimarySequenceSchema).default([]),
});

export type ClusterModule = z.output<typeof ClusterModuleSchema>;
export type ClusterReadout = z.output<typeof ClusterReadoutSchema>;
export type ClusterPayload = z.output<typeof ClusterPayloadSchema>;
export type ClusterPrimarySequence = z.output<typeof ClusterPrimarySequenceSchema>;
