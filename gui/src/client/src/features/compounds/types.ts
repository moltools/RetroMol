import { z } from "zod";

export const CompoundOptionSchema = z.object({
  name: z.string(),
  smiles: z.string(),
  databaseName: z.string(),
  // Backend sends this as a DB row id (number); normalize to string for the UI.
  databaseIdentifier: z.union([z.string(), z.number()]).transform(String),
});
export type CompoundOption = z.output<typeof CompoundOptionSchema>;

export const SearchCompoundRespSchema = z.object({
  rows: z.array(CompoundOptionSchema).default([]),
  rowCount: z.number().int().nonnegative().optional(),
});
