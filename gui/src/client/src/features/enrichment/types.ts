import { z } from "zod";

export const EntryTypeSchema = z.enum(["compound", "bgc"]);
export type EntryType = z.output<typeof EntryTypeSchema>;

export const EntrySourceSchema = z.object({
  name: z.string(),
  databaseName: z.string(),
  url: z.string().nullable(),
});

export const SearchEntrySchema = z.object({
  id: z.string(),
  name: z.string(),
  url: z.string().nullable(),
  type: EntryTypeSchema,
  sources: z.array(EntrySourceSchema),
});
export type SearchEntry = z.output<typeof SearchEntrySchema>;

export const EntrySearchRespSchema = z.object({
  results: z.array(SearchEntrySchema),
});
export type EntrySearchResp = z.output<typeof EntrySearchRespSchema>;

export const EnrichmentResultSchema = z.object({
  termId: z.string(),
  category: z.string().nullable(),
  rank: z.string().nullable(),
  label: z.string(),
  selectedWithTerm: z.number().int().nonnegative(),
  selectedTotal: z.number().int().nonnegative(),
  backgroundWithTerm: z.number().int().nonnegative(),
  backgroundTotal: z.number().int().nonnegative(),
  foldEnrichment: z.number().nullable(),
  direction: z.enum(["enriched", "depleted"]),
  pValue: z.number(),
  qValue: z.number(),
});
export type EnrichmentResult = z.output<typeof EnrichmentResultSchema>;

export const EnrichmentAnalysisRespSchema = z.object({
  results: z.array(EnrichmentResultSchema),
});
export type EnrichmentAnalysisResp = z.output<typeof EnrichmentAnalysisRespSchema>;
