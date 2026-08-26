import { z } from "zod";

export const CountSchema = z.object({
  label: z.string(),
  count: z.number().int().nonnegative(),
});
export type Count = z.output<typeof CountSchema>;

export const DatabaseStatsRespSchema = z.object({
  totalEntries: z.number().int().nonnegative(),
  countsByType: z.array(CountSchema),
  sequenceLengthMin: z.number().int().nonnegative(),
  sequenceLengthMax: z.number().int().nonnegative(),
  sequenceLengthAvg: z.number().nonnegative(),
  uniqueBlockCount: z.number().int().nonnegative(),
  withSourceUrlCount: z.number().int().nonnegative(),
  withoutSourceUrlCount: z.number().int().nonnegative(),
});
export type DatabaseStatsResp = z.output<typeof DatabaseStatsRespSchema>;

export const AnnotationStatsRespSchema = z.object({
  withAnnotationCount: z.number().int().nonnegative(),
  withoutAnnotationCount: z.number().int().nonnegative(),
  countsByCategory: z.array(CountSchema),
  phylogenyTypeCounts: z.array(CountSchema),
  phylogenyGenusCounts: z.array(CountSchema),
  phylogenySpeciesCounts: z.array(CountSchema),
  biosyntheticClassCounts: z.array(CountSchema),
  chemicalClassPathwayCounts: z.array(CountSchema),
  chemicalClassSuperclassCounts: z.array(CountSchema),
  chemicalClassClassCounts: z.array(CountSchema),
  bioactivityAtcCounts: z.array(CountSchema),
  bioactivityMaxPhaseCounts: z.array(CountSchema),
  bioactivityBiologicalRoleCounts: z.array(CountSchema),
  bioactivityChemicalRoleCounts: z.array(CountSchema),
});
export type AnnotationStatsResp = z.output<typeof AnnotationStatsRespSchema>;

export const EntryAnnotationSchema = z.object({
  id: z.string(),
  category: z.string(),
  rank: z.string().nullable(),
  label: z.string(),
  externalId: z.string().nullable(),
  url: z.string().nullable(),
});
export type EntryAnnotation = z.output<typeof EntryAnnotationSchema>;

export const EntryAnnotationsRespSchema = z.object({
  results: z.array(EntryAnnotationSchema),
});
export type EntryAnnotationsResp = z.output<typeof EntryAnnotationsRespSchema>;
