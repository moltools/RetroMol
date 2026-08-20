import { getJson } from "../http";
import { AnnotationStatsRespSchema, DatabaseStatsRespSchema, type AnnotationStatsResp, type DatabaseStatsResp } from "./types";

export async function getDatabaseStats(signal?: AbortSignal): Promise<DatabaseStatsResp> {
  return getJson("/api/databaseStats", DatabaseStatsRespSchema, signal);
}

export async function getAnnotationStats(signal?: AbortSignal): Promise<AnnotationStatsResp> {
  return getJson("/api/annotationStats", AnnotationStatsRespSchema, signal);
}
