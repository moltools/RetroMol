import { getJson, postJson } from "../http";
import {
  SubmitDiscoveryQueryRespSchema,
  GetDiscoveryQueryResultRespSchema,
  type DiscoveryQueryItem,
  type DiscoveryQueryPayload,
} from "../session/types";
import {
  MonomerNameSearchRespSchema,
  DiscoveryQueryReqSchema,
  DiscoveryQueryRespSchema,
  DiscoveryMsaReqSchema,
  DiscoveryMsaRespSchema,
  DiscoveryCompareTanimotoReqSchema,
  DiscoveryCompareTanimotoRespSchema,
  SubmitDiscoveryQueryReqSchema,
  type DiscoveryQueryReq,
  type MonomerNameOption,
  type DiscoveryQueryResp,
  type DiscoveryMsaReq,
  type DiscoveryMsaResp,
  type DiscoveryCompareTanimotoReq,
  type DiscoveryCompareTanimotoResp,
  type SubmitDiscoveryQueryReq,
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

export async function runDiscoveryMsa(payload: DiscoveryMsaReq, signal?: AbortSignal): Promise<DiscoveryMsaResp> {
  const validated = DiscoveryMsaReqSchema.parse(payload);
  return postJson("/api/discoveryMsa", validated, DiscoveryMsaRespSchema, signal);
}

export async function runDiscoveryCompareTanimoto(
  payload: DiscoveryCompareTanimotoReq,
  signal?: AbortSignal
): Promise<DiscoveryCompareTanimotoResp> {
  const validated = DiscoveryCompareTanimotoReqSchema.parse(payload);
  return postJson("/api/discoveryCompareTanimoto", validated, DiscoveryCompareTanimotoRespSchema, signal);
}

// Fire-and-forget: the created item (status "queued") comes back immediately, and its
// status/payload updates flow through the existing SSE + getSession mechanism from
// there -- same pattern as jobs/api.ts's submitCompoundJob, just returning the item
// since (unlike a compound) it isn't created client-side first.
export async function submitDiscoveryQuery(payload: SubmitDiscoveryQueryReq): Promise<DiscoveryQueryItem> {
  const validated = SubmitDiscoveryQueryReqSchema.parse(payload);
  const data = await postJson("/api/submitDiscoveryQuery", validated, SubmitDiscoveryQueryRespSchema);
  return data.item;
}

// payload is stripped from session items by getSession -- this fetches a discovery
// query's actual computed payload on demand, the same way reconstructCompound/
// getClusterReadout do for the other item kinds.
export async function getDiscoveryQueryResult(
  sessionId: string,
  itemId: string,
  signal?: AbortSignal
): Promise<DiscoveryQueryPayload | null> {
  const data = await postJson("/api/getDiscoveryQueryResult", { sessionId, itemId }, GetDiscoveryQueryResultRespSchema, signal);
  return data.payload;
}
