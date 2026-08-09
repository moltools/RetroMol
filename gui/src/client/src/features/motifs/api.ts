import { getJson } from "../http";
import { MotifStructuresRespSchema, type MotifStructures } from "./types";

// The whole name -> SMILES vocabulary, fetched once and cached by the caller
// (see MotifHoverCard) -- it's small and effectively static for the life of
// the server process, same rationale as searchMonomerNames' rule list.
export async function fetchMotifStructures(signal?: AbortSignal): Promise<MotifStructures> {
  const data = await getJson("/api/motifStructures", MotifStructuresRespSchema, signal);
  return data.structures;
}
