import { postJson } from "../http";
import { ClusterPayload, ClusterPrimarySequence, GetClusterReadoutRespSchema, ReconstructGeneClusterRespSchema } from "./types";

// A gene cluster's parsed linear module readout(s) are computed once at submit
// time and stashed server-side on the item -- but `payload` is deliberately
// stripped from items before they reach the client as part of the session (see
// strip_property_from_dict in session_store.py), same as a compound's `payload`
// only ever reaching the client via /api/reconstructCompound. This is the
// gene-cluster equivalent of that endpoint.
export async function getClusterReadout(
  sessionId: string,
  itemId: string,
  signal?: AbortSignal
): Promise<ClusterPayload> {
  const resp = await postJson(
    "/api/getClusterReadout",
    { sessionId, itemId },
    GetClusterReadoutRespSchema,
    signal
  );
  return resp.data;
}

// Maps a gene cluster's readout(s) onto the same matching-rule vocabulary a
// compound's primary sequence is drawn from -- the gene-cluster equivalent of
// reconstructCompound, so a BGC upload can seed a discovery query the same way.
export async function reconstructGeneCluster(
  sessionId: string,
  itemId: string,
  signal?: AbortSignal
): Promise<ClusterPrimarySequence[]> {
  const resp = await postJson(
    "/api/reconstructGeneCluster",
    { sessionId, itemId },
    ReconstructGeneClusterRespSchema,
    signal
  );
  return resp.data;
}
