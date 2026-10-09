import { getJson } from "../http";
import { SearchCompoundRespSchema, CompoundOption } from "./types";

export async function searchCompounds(
  query: string,
  limit: number = 10,
  signal?: AbortSignal
): Promise<CompoundOption[]> {
  const params = new URLSearchParams({ q: query, limit: String(limit) });
  const data = await getJson(`/api/searchCompound?${params.toString()}`, SearchCompoundRespSchema, signal);
  return data.rows;
}
