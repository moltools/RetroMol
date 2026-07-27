import { getJson, postJson } from "../http";
import {
  MonomerNameSearchRespSchema,
  DiscoveryQueryReqSchema,
  DiscoveryQueryRespSchema,
  type DiscoveryQueryReq,
  type MonomerNameOption,
  type DiscoveryQueryResp,
} from "./types";

export async function searchMonomerNames(
  query: string,
  limit: number = 10,
  signal?: AbortSignal
): Promise<MonomerNameOption[]> {
  const params = new URLSearchParams({ q: query, limit: String(limit) });
  const data = await getJson(`/api/discoveryMonomerNames?${params.toString()}`, MonomerNameSearchRespSchema, signal);
  return data.rows;
}

export async function runDiscoveryQuery(
  payload: DiscoveryQueryReq,
  signal?: AbortSignal
): Promise<DiscoveryQueryResp> {
  const validated = DiscoveryQueryReqSchema.parse(payload);
  return postJson("/api/discoveryQuery", validated, DiscoveryQueryRespSchema, signal);
}
