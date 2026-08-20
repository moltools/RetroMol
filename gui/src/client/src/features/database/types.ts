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
  topGenera: z.array(CountSchema),
  chemicalClassCounts: z.array(CountSchema),
});
export type AnnotationStatsResp = z.output<typeof AnnotationStatsRespSchema>;
