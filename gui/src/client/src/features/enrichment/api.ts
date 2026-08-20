import { getJson, postJson } from "../http";
import {
  EntrySearchRespSchema,
  EnrichmentAnalysisRespSchema,
  type EntrySearchResp,
  type EnrichmentAnalysisResp,
  type EntryType,
} from "./types";

export const MAX_ENTRY_SEARCH_RESULTS = 100;
export const MAX_ENRICHMENT_SELECTION = 100;

export async function searchEntries(
  query: string,
  entryType: EntryType | "all" = "all",
  signal?: AbortSignal
): Promise<EntrySearchResp> {
  const params = new URLSearchParams({ q: query, type: entryType, limit: String(MAX_ENTRY_SEARCH_RESULTS) });
  return getJson(`/api/entrySearch?${params.toString()}`, EntrySearchRespSchema, signal);
}

export async function runEnrichmentAnalysis(
  entryIds: string[],
  signal?: AbortSignal
): Promise<EnrichmentAnalysisResp> {
  return postJson("/api/enrichmentAnalysis", { entryIds }, EnrichmentAnalysisRespSchema, signal);
}
