import { getJson } from "../http";
import { MotifStructuresRespSchema, type MotifStructuresData } from "./types";

// The whole name -> SMILES vocabulary (plus which names are ambiguous -- see
// MotifStructuresData), fetched once and cached by the caller (see
// MotifHoverCard) -- it's small and effectively static for the life of the
// server process, same rationale as searchMonomerNames' rule list.
export async function fetchMotifStructures(signal?: AbortSignal): Promise<MotifStructuresData> {
  return getJson("/api/motifStructures", MotifStructuresRespSchema, signal);
}
