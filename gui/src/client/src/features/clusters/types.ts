import { z } from "zod";
import { PrimarySequenceItemSchema } from "../reconstruction/types";

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
  // Resolved via pmp.yml's "pks.extender"/"pks.kr_stereochemistry" axes (see
  // retromol_antismash/modules.py's PKSAnatomy) -- null when unresolved (e.g.
  // trans-AT modules, or an unmapped antiSMASH substrate call). These aren't
  // enough on their own to know the module's *specific* token (e.g. "A2",
  // "B^R2") -- that also depends on which name actually exists in mxn.yml, so
  // prefer ReconstructGeneClusterResp's primary_sequence for that instead of
  // recomputing it client-side from these flags.
  extender_digit: z.string().nullable().optional(),
  beta_stereo: z.string().nullable().optional(),
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
// retromol_antismash.modules.bgc_primary_sequence on the backend).
//
// `modules[i]` is the module `primary_sequence[i]`/`tokens[i]` came from -- all
// three arrays are index-aligned, in biosynthetic order (NOT the same order as
// ClusterReadoutSchema.modules, which is raw insertion order: see the comment on
// gui/src/server/routes/jobs.py's run_gene_cluster_reconstruction).
//
// Each primary_sequence entry is a [name, [moduleIndex]] pair -- the same tuple
// shape as a compound's PrimarySequenceItem, so the same click-to-highlight chip
// components (PrimarySequenceChips, SequenceEditor) work unmodified. "tags" here
// is always a single-element array holding the entry's index into `modules`
// (not an atom index, unlike the compound case) -- use it to link a readout chip
// back to the gene/domains it came from.
export const ClusterPrimarySequenceSchema = z.object({
  id: z.string(),
  primary_sequence: z.array(PrimarySequenceItemSchema),
  tokens: z.array(z.array(z.string())).default([]),
  modules: z.array(ClusterModuleSchema).default([]),
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
