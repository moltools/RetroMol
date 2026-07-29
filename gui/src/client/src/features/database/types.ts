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
  sequenceLengthBuckets: z.array(CountSchema),
  topBlocks: z.array(CountSchema),
  uniqueBlockCount: z.number().int().nonnegative(),
  fullyResolvedCount: z.number().int().nonnegative(),
  hasUnknownBlockCount: z.number().int().nonnegative(),
  withSourceUrlCount: z.number().int().nonnegative(),
  withoutSourceUrlCount: z.number().int().nonnegative(),
});
export type DatabaseStatsResp = z.output<typeof DatabaseStatsRespSchema>;
