import { getJson } from "../http";
import { DatabaseStatsRespSchema, type DatabaseStatsResp } from "./types";

export async function getDatabaseStats(signal?: AbortSignal): Promise<DatabaseStatsResp> {
  return getJson("/api/databaseStats", DatabaseStatsRespSchema, signal);
}
