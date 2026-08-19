import { z } from "zod";

export const PrimarySequenceItemSchema = z.tuple([z.string(), z.array(z.number())]);
export type PrimarySequenceItem = z.output<typeof PrimarySequenceItemSchema>;

export const ReconstructionSchema = z.object({
  tagged_input_smiles: z.string(),
  // Null when RetroMol's hardcoded fusion chemistry couldn't rebuild a backbone
  // for this path -- see `backbone_warning` in that case. The primary sequence is
  // still meaningful (it's derived directly from the source structure), so a
  // candidate is kept even without a backbone.
  tagged_backbone_smiles: z.string().nullable(),
  primary_sequence: z.array(PrimarySequenceItemSchema),
  backbone_warning: z.string().nullable().default(null),
  // False when `primary_sequence` isn't a genuine biosynthetic order -- RetroMol
  // identified every building block individually but couldn't connect them into a
  // single path (e.g. a branched or cyclic assembly). Render as an unordered set,
  // not a sequence.
  ordered: z.boolean().default(true),
  // Names of this compound's tailoring events (glycosylation, methylation, ...) --
  // real identified content that isn't part of any chain, so it's never displayed
  // as part of primary_sequence, but should still count toward a discovery query's
  // fingerprint (see SubmitDiscoveryQueryReqSchema.extraFingerprintTokens). Always
  // empty for `data` (reconstruct_linear_readout candidates don't carry this);
  // only ever populated on `dbMatchingSequences` entries.
  extra_fingerprint_tokens: z.array(z.string()).default([]),
});
export type Reconstruction = z.output<typeof ReconstructionSchema>;

export const ReconstructCompoundRespSchema = z.object({
  ok: z.boolean().optional(),
  status: z.string().optional(),
  data: z.array(ReconstructionSchema).default([]),
  // Candidates read directly off result.linear_readout.paths -- the same
  // representation the database is actually populated with, unlike `data` above
  // (reconstruct_linear_readout's candidates, which apply backbone-reconstruction
  // eligibility/orientation filtering that can diverge from what's stored). Query
  // with one of these, not `data`, to get a fingerprint comparable to the database.
  dbMatchingSequences: z.array(ReconstructionSchema).default([]),
});
