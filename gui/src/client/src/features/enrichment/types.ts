import { z } from "zod";

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

// q-value threshold below which a term is flagged significant in the results table.
export const Q_VALUE_SIGNIFICANT = 0.05;
