import { getJson } from "../http";
import {
  AnnotationStatsRespSchema,
  DatabaseStatsRespSchema,
  EntryAnnotationsRespSchema,
  type AnnotationStatsResp,
  type DatabaseStatsResp,
  type EntryAnnotationsResp,
} from "./types";

export async function getDatabaseStats(signal?: AbortSignal): Promise<DatabaseStatsResp> {
  return getJson("/api/databaseStats", DatabaseStatsRespSchema, signal);
}

export async function getAnnotationStats(signal?: AbortSignal): Promise<AnnotationStatsResp> {
  return getJson("/api/annotationStats", AnnotationStatsRespSchema, signal);
}

export async function getEntryAnnotations(entryId: string, signal?: AbortSignal): Promise<EntryAnnotationsResp> {
  const params = new URLSearchParams({ entryId });
  return getJson(`/api/entryAnnotations?${params.toString()}`, EntryAnnotationsRespSchema, signal);
}
