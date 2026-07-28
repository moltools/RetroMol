import { z } from "zod";

export const MonomerNameOptionSchema = z.object({
  name: z.string(),
});
export type MonomerNameOption = z.output<typeof MonomerNameOptionSchema>;

export const MonomerNameSearchRespSchema = z.object({
  rows: z.array(MonomerNameOptionSchema).default([]),
  rowCount: z.number().int().nonnegative().optional(),
});

export const DiscoveryEntryTypeSchema = z.enum(["compound", "bgc", "both"]);
export type DiscoveryEntryType = z.output<typeof DiscoveryEntryTypeSchema>;

// What the normalized alignment score is computed relative to. "subsequence"
// (default) only considers the query's own length, so a short query matching
// well as a subsequence of a much longer candidate can still score near 100% —
// intentional, since that mode is meant to surface subsequence matches.
// "longest_sequence" instead normalizes against whichever of the two sequences
// is longer, penalizing that case and favoring matches between similarly-sized
// sequences. Add new modes here as they're added on the backend.
export const DiscoveryScoreModeSchema = z.enum(["subsequence", "longest_sequence"]);
export type DiscoveryScoreMode = z.output<typeof DiscoveryScoreModeSchema>;

export const DISCOVERY_SCORE_MODE_OPTIONS: { value: DiscoveryScoreMode; label: string }[] = [
  { value: "subsequence", label: "Favor subsequence matches" },
  { value: "longest_sequence", label: "Favor equally-sized matches" },
];

export const DiscoveryQueryReqSchema = z.object({
  primarySequence: z.array(z.string().min(1)).min(1, "Sequence must have at least one block"),
  entryType: DiscoveryEntryTypeSchema,
  scoreMode: DiscoveryScoreModeSchema,
  n: z.number().int().min(1).max(1000),
  topX: z.number().int().min(1),
});
export type DiscoveryQueryReq = z.output<typeof DiscoveryQueryReqSchema>;

export const DiscoveryResultSchema = z.object({
  entryId: z.string(),
  name: z.string(),
  url: z.string().nullable(),
  type: z.enum(["compound", "bgc"]),
  fingerprintSimilarity: z.number(),
  primarySequence: z.array(z.string()),
  inverted: z.boolean(),
  alignmentScore: z.number(),
  normalizedAlignmentScorePct: z.number(),
  alignedQuery: z.array(z.string().nullable()),
  alignedTarget: z.array(z.string().nullable()),
  // Per-column Tanimoto similarity (0-1) between alignedQuery[i] and alignedTarget[i],
  // from the same scoring matrix used to compute the alignment; null wherever either
  // side is a gap. Index-aligned with alignedQuery/alignedTarget.
  alignedSimilarity: z.array(z.number().nullable()),
});
export type DiscoveryResult = z.output<typeof DiscoveryResultSchema>;

export const DiscoveryQueryRespSchema = z.object({
  ok: z.boolean().optional(),
  scoreMode: DiscoveryScoreModeSchema.optional(),
  querySequence: z.array(z.string()),
  selfScore: z.number(),
  candidatesConsidered: z.number().int().nonnegative(),
  candidatesSkipped: z.number().int().nonnegative(),
  results: z.array(DiscoveryResultSchema),
});
export type DiscoveryQueryResp = z.output<typeof DiscoveryQueryRespSchema>;

export const DiscoveryMsaReqSchema = z.object({
  querySequence: z.array(z.string().min(1)).min(1),
  sequences: z
    .array(
      z.object({
        id: z.string().min(1),
        sequence: z.array(z.string().min(1)).min(1),
      })
    )
    .min(1),
});
export type DiscoveryMsaReq = z.output<typeof DiscoveryMsaReqSchema>;

// "query" is a reserved id identifying the query's own row in the response.
export const DiscoveryMsaRowSchema = z.object({
  id: z.string(),
  alignedSequence: z.array(z.string().nullable()),
  // Per-column Tanimoto similarity (0-1) against the query's own row at that same MSA
  // column; all null for the query's row itself. Index-aligned with alignedSequence.
  similarityToQuery: z.array(z.number().nullable()),
});
export type DiscoveryMsaRow = z.output<typeof DiscoveryMsaRowSchema>;

export const DiscoveryMsaRespSchema = z.object({
  ok: z.boolean().optional(),
  rows: z.array(DiscoveryMsaRowSchema),
});
export type DiscoveryMsaResp = z.output<typeof DiscoveryMsaRespSchema>;
