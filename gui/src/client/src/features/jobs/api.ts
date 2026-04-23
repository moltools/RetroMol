import { postJson } from "../http";
import type { WorkspaceImportDeps, NewCompoundJob } from "./types";
import type { Session, SessionItem, CompoundItem, ClusterItem } from "../session/types";
import { saveSession } from "../session/api";
import { z } from "zod";

export const MAX_ITEMS = 50;

const SubmitJobRespSchema = z.object({
  ok: z.boolean(),
  elapsed_ms: z.number().int().nonnegative(),
  status: z.string().optional(),
}).partial();

export async function submitCompoundJob(
  sessionId: string,
  item: CompoundItem,
): Promise<void> {
  await postJson(
    "/api/submitCompound",
    {
      sessionId,
      itemId: item.id,
      name: item.name,
      smiles: item.smiles,
      matchStereochemistry: item.matchStereochemistry,
    },
    SubmitJobRespSchema
  );
};

export async function submitClusterJob(
  sessionId: string,
  item: ClusterItem,
): Promise<void> {
  await postJson(
    "/api/submitCluster",
    {
      sessionId,
      itemId: item.id,
      name: item.name,
      fileContent: item.fileContent,
    },
    SubmitJobRespSchema
  );
};

// Batch compound import
export async function importCompoundsBatch(
  deps: WorkspaceImportDeps,
  compounds: NewCompoundJob[],
): Promise<SessionItem[]> {
  const { pushNotification, setSession, sessionId } = deps;

  if (!compounds.length) {
    pushNotification("No compounds to import", "warning");
    return [];
  };

  let nextSession: Session | null = null;
  let newItems: SessionItem[] = [];

  // Update local session (queued items)
  setSession((prev) => {
    const existingCount = prev.items.length;
    const remainingSlots = MAX_ITEMS - existingCount;

    if (remainingSlots <= 0) {
      pushNotification(`Session already has maximum of ${MAX_ITEMS} items`, "warning");
      nextSession = prev;
      newItems = [];
      return prev;
    };

    const limited = compounds.length > remainingSlots ? compounds.slice(0, remainingSlots) : compounds;

    if (limited.length < compounds.length) {
      pushNotification(`Only importing ${limited.length} compounds to avoid exceeding maximum of ${MAX_ITEMS} items`, "warning");
    };

    const createdItems: SessionItem[] = limited.map(({ name, smiles, matchStereochemistry }) => ({
      id: crypto.randomUUID(),
      kind: "compound",
      name,
      smiles,
      matchStereochemistry,
      status: "queued",
      errorMessage: null,
      updatedAt: Date.now(),
      // optional fields
      score: null,
    }));

    const updated: Session = { ...prev, items: [...prev.items, ...createdItems] };

    nextSession = updated;
    newItems = createdItems;
    return updated;
  });

  if (!nextSession || newItems.length === 0) return [];

  // Persist session BEFORE submitting jobs
  try {
    await saveSession(nextSession);
  } catch (err) {
    const msg = err instanceof Error ? err.message : String(err);
    pushNotification(`Failed to save session before importing compounds: ${msg}`, "error");

    const newIds = new Set(newItems.map((it) => it.id));

    setSession((prev) => ({
      ...prev,
      items: prev.items.map((it) =>
        newIds.has(it.id)
          ? {
              ...it,
              status: "error",
              errorMessage: "Failed to save session before importing compound",
              updatedAt: Date.now(),
            }
          : it
      )
    }));

    return [];
  };

  // Submit jobs sequentially
  for (const item of newItems) {
    try {
      await submitCompoundJob(sessionId, item as CompoundItem);
    } catch (err) {
      const msg = err instanceof Error ? err.message : String(err);
      pushNotification(`Failed to submit job for compound "${item.name}": ${msg}`, "error");

      // Mark item as error
      setSession((prev) => ({
        ...prev,
        items: prev.items.map((it) =>
          it.id === item.id
            ? {
                ...it,
                status: "error",
                errorMessage: `Failed to submit job: ${msg}`,
                updatedAt: Date.now(),
              }
            : it
        )
      }));
    };
  };

  return newItems;
};

// Single compound import wrapper
export async function importCompound(
  deps: WorkspaceImportDeps,
  payload: NewCompoundJob,
): Promise<SessionItem | null> {
  const items = await importCompoundsBatch(deps, [payload]);
  return items[0] ?? null;
};

// Batch cluster import
export async function importClustersBatch(
  deps: WorkspaceImportDeps,
  clusters: { name: string; fileContent: string }[],
): Promise<SessionItem[]> {
  const { pushNotification, setSession, sessionId } = deps;

  if (!clusters.length) {
    pushNotification("No clusters to import", "warning");
    return [];
  };

  let nextSession: Session | null = null;
  let newItems: SessionItem[] = [];

  // Update local session (queued items)
  setSession((prev) => {
    const existingCount = prev.items.length;
    const remainingSlots = MAX_ITEMS - existingCount;

    if (remainingSlots <= 0) {
      pushNotification(`Session already has maximum of ${MAX_ITEMS} items`, "warning");
      nextSession = prev;
      newItems = [];
      return prev;
    };

    const limited = clusters.length > remainingSlots ? clusters.slice(0, remainingSlots) : clusters;

    if (limited.length < clusters.length) {
      pushNotification(`Only importing ${limited.length} clusters to avoid exceeding maximum of ${MAX_ITEMS} items`, "warning");
    };

    const createdItems: SessionItem[] = limited.map(({ name, fileContent }) => ({
      id: crypto.randomUUID(),
      kind: "cluster",
      name,
      fileContent,
      status: "queued",
      errorMessage: null,
      updatedAt: Date.now(),
    }));

    const updated: Session = { ...prev, items: [...prev.items, ...createdItems] };

    nextSession = updated;
    newItems = createdItems;
    return updated;
  });

  if (!nextSession || newItems.length === 0) return [];

  // Persist session BEFORE submitting jobs
  try {
    await saveSession(nextSession);
  } catch (err) {
    const msg = err instanceof Error ? err.message : String(err);
    pushNotification(`Failed to save session before importing clusters: ${msg}`, "error");

    const newIds = new Set(newItems.map((it) => it.id));

    setSession((prev) => ({
      ...prev,
      items: prev.items.map((it) =>
        newIds.has(it.id)
          ? {
              ...it,
              status: "error",
              errorMessage: "Failed to save session before importing cluster",
              updatedAt: Date.now(),
            }
          : it
      )
    }));

    return [];
  };

  // Submit jobs sequentially
  for (const item of newItems) {
    try {
      await submitClusterJob(sessionId, item as ClusterItem);
    } catch (err) {
      const msg = err instanceof Error ? err.message : String(err);
      pushNotification(`Failed to submit job for cluster "${item.name}": ${msg}`, "error");

      // Mark item as error
      setSession((prev) => ({
        ...prev,
        items: prev.items.map((it) =>
          it.id === item.id
            ? {
                ...it,
                status: "error",
                errorMessage: `Failed to submit job: ${msg}`,
                updatedAt: Date.now(),
              }
            : it
        )
      }));
    };
  };

  return newItems;
};