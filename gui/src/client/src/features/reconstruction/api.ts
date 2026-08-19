import { postJson } from "../http";
import { saveSession } from "../session/api";
import type { Session } from "../session/types";
import { ReconstructCompoundRespSchema, Reconstruction, PrimarySequenceItem } from "./types";

export interface ReconstructCompoundResult {
  reconstructions: Reconstruction[];
  // Candidates built to match the database's own recipe -- prefer these for a
  // discovery query. See ReconstructCompoundRespSchema.dbMatchingSequences.
  dbMatchingSequences: Reconstruction[];
}

export async function reconstructCompound(
  sessionId: string,
  itemId: string,
  signal?: AbortSignal
): Promise<ReconstructCompoundResult> {
  const data = await postJson(
    "/api/reconstructCompound",
    { sessionId, itemId },
    ReconstructCompoundRespSchema,
    signal
  );
  return { reconstructions: data.data, dbMatchingSequences: data.dbMatchingSequences };
}

function withEditedPrimarySequences(
  session: Session,
  itemId: string,
  mutate: (current: Record<string, PrimarySequenceItem[]>) => Record<string, PrimarySequenceItem[]>,
): Session {
  return {
    ...session,
    items: session.items.map((it) => {
      if (it.id !== itemId || it.kind !== "compound") return it;
      return { ...it, editedPrimarySequences: mutate(it.editedPrimarySequences ?? {}) };
    }),
  };
}

// Persists user-edited primary sequence(s) for one compound item, keyed by
// reconstruction index. Only the given indices are touched -- any other saved
// overrides on this item (or edits on other items) are carried over as-is.
export async function saveEditedPrimarySequences(
  session: Session,
  itemId: string,
  overrides: Record<number, PrimarySequenceItem[]>,
): Promise<Session> {
  const nextSession = withEditedPrimarySequences(session, itemId, (current) => {
    const next = { ...current };
    for (const [idx, sequence] of Object.entries(overrides)) {
      next[idx] = sequence;
    }
    return next;
  });
  await saveSession(nextSession);
  return nextSession;
}

// Clears a saved edit for one reconstruction index. We never persist a copy of
// the algorithm's own output, so "reverting" is just deleting the override --
// the next render falls back to the (always freshly-fetched) parsed sequence.
export async function revertEditedPrimarySequence(
  session: Session,
  itemId: string,
  reconstructionIndex: number,
): Promise<Session> {
  const nextSession = withEditedPrimarySequences(session, itemId, (current) => {
    const { [String(reconstructionIndex)]: _removed, ...rest } = current;
    return rest;
  });
  await saveSession(nextSession);
  return nextSession;
}
