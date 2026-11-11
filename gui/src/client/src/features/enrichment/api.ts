import { postJson } from "../http";
import { EnrichmentAnalysisRespSchema, type EnrichmentAnalysisResp } from "./types";

export const MAX_ENRICHMENT_SELECTION = 100;

export async function runEnrichmentAnalysis(
  entryIds: string[],
  signal?: AbortSignal
): Promise<EnrichmentAnalysisResp> {
  return postJson("/api/enrichmentAnalysis", { entryIds }, EnrichmentAnalysisRespSchema, signal);
}
