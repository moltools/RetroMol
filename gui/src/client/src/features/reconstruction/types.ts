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
});
export type Reconstruction = z.output<typeof ReconstructionSchema>;

export const ReconstructCompoundRespSchema = z.object({
  ok: z.boolean().optional(),
  status: z.string().optional(),
  data: z.array(ReconstructionSchema).default([]),
});
