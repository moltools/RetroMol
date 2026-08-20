import { z } from "zod";

// name -> SMILES to depict for that motif. Absent keys mean "no matching rule
// for this name" (unidentified "X" blocks, hand-edited names, PK_GROUP_TOKENS).
export const MotifStructuresRespSchema = z.object({
  structures: z.record(z.string()),
  // Names a matching rule doesn't uniquely identify -- e.g. "glycosylation" alone
  // names dozens of rules, one per distinct sugar. `structures[name]` for one of
  // these is just whichever rule happened to be registered first, not necessarily
  // the one that actually matched a given occurrence of that name.
  ambiguousNames: z.array(z.string()).default([]),
});

export type MotifStructures = z.output<typeof MotifStructuresRespSchema>["structures"];
export type MotifStructuresData = z.output<typeof MotifStructuresRespSchema>;
