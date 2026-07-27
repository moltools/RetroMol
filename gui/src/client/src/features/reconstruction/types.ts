import { z } from "zod";

export const PrimarySequenceItemSchema = z.tuple([z.string(), z.array(z.number())]);
export type PrimarySequenceItem = z.output<typeof PrimarySequenceItemSchema>;

export const ReconstructionSchema = z.object({
  tagged_input_smiles: z.string(),
  tagged_backbone_smiles: z.string(),
  primary_sequence: z.array(PrimarySequenceItemSchema),
});
export type Reconstruction = z.output<typeof ReconstructionSchema>;

export const ReconstructCompoundRespSchema = z.object({
  ok: z.boolean().optional(),
  status: z.string().optional(),
  data: z.array(ReconstructionSchema).default([]),
});
