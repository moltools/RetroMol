import { getJson } from "../http";
import { AnnotationTermsRespSchema, BrowseEntriesRespSchema } from "./types";
import type { EntryType } from "../enrichment/types";

function buildParams(entryType: EntryType | "all", termId: string | null): URLSearchParams {
  const params = new URLSearchParams();
  if (entryType !== "all") params.set("type", entryType);
  if (termId) params.set("termId", termId);
  return params;
}

export async function getAnnotationTerms(category?: string, signal?: AbortSignal) {
  const params = new URLSearchParams();
  if (category) params.set("category", category);
  const qs = params.toString();
  return getJson(`/api/annotationTerms${qs ? `?${qs}` : ""}`, AnnotationTermsRespSchema, signal);
}

export async function getBrowseEntries(
  entryType: EntryType | "all" = "all",
  termId: string | null = null,
  signal?: AbortSignal
) {
  const params = buildParams(entryType, termId);
  return getJson(`/api/browseEntries?${params.toString()}`, BrowseEntriesRespSchema, signal);
}

export function browseEntriesExportUrl(entryType: EntryType | "all", termId: string | null): string {
  const params = buildParams(entryType, termId);
  return `/api/browseEntries.tsv?${params.toString()}`;
}
