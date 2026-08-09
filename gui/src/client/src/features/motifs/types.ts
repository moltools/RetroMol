import { z } from "zod";

// name -> SMILES to depict for that motif. Absent keys mean "no matching rule
// for this name" (unidentified "X" blocks, hand-edited names, PK_GROUP_TOKENS).
export const MotifStructuresRespSchema = z.object({
  structures: z.record(z.string()),
});

export type MotifStructures = z.output<typeof MotifStructuresRespSchema>["structures"];
