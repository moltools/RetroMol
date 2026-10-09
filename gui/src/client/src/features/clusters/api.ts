import { postJson } from "../http";
import { saveSession } from "../session/api";
import type { Session } from "../session/types";
import type { PrimarySequenceItem } from "../reconstruction/types";
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

function withEditedClusterPrimarySequences(
  session: Session,
  itemId: string,
  mutate: (current: Record<string, PrimarySequenceItem[]>) => Record<string, PrimarySequenceItem[]>,
): Session {
  return {
    ...session,
    items: session.items.map((it) => {
      if (it.id !== itemId || it.kind !== "cluster") return it;
      return { ...it, editedPrimarySequences: mutate(it.editedPrimarySequences ?? {}) };
    }),
  };
}

// Persists a user-edited primary sequence for one gene-cluster region, keyed by
// its index in the item's readouts -- mirrors reconstruction/api.ts's
// saveEditedPrimarySequences for compounds.
export async function saveEditedClusterPrimarySequence(
  session: Session,
  itemId: string,
  regionIndex: number,
  sequence: PrimarySequenceItem[],
): Promise<Session> {
  const nextSession = withEditedClusterPrimarySequences(session, itemId, (current) => ({
    ...current,
    [String(regionIndex)]: sequence,
  }));
  await saveSession(nextSession);
  return nextSession;
}

// Clears a saved edit for one region -- mirrors revertEditedPrimarySequence.
export async function revertEditedClusterPrimarySequence(
  session: Session,
  itemId: string,
  regionIndex: number,
): Promise<Session> {
  const nextSession = withEditedClusterPrimarySequences(session, itemId, (current) => {
    const { [String(regionIndex)]: _removed, ...rest } = current;
    return rest;
  });
  await saveSession(nextSession);
  return nextSession;
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
