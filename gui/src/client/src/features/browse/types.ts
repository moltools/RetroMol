import { z } from "zod";
import { EntryTypeSchema, EntrySourceSchema } from "../enrichment/types";

export const AnnotationTermSchema = z.object({
  id: z.string(),
  category: z.string(),
  rank: z.string().nullable(),
  label: z.string(),
  parentId: z.string().nullable(),
});
export type AnnotationTerm = z.output<typeof AnnotationTermSchema>;

export const AnnotationTermsRespSchema = z.object({
  terms: z.array(AnnotationTermSchema),
});

export const BrowseEntrySchema = z.object({
  id: z.string(),
  type: EntryTypeSchema,
  name: z.string(),
  url: z.string().nullable(),
  smiles: z.string().nullable(),
  sources: z.array(EntrySourceSchema),
  phylogenyType: z.string().nullable(),
  genus: z.string().nullable(),
  species: z.string().nullable(),
  chemicalClasses: z.array(z.string()),
});
export type BrowseEntry = z.output<typeof BrowseEntrySchema>;

export const BrowseEntriesRespSchema = z.object({
  entries: z.array(BrowseEntrySchema),
});
export type BrowseEntriesResp = z.output<typeof BrowseEntriesRespSchema>;
