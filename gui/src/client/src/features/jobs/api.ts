import { postJson } from "../http";
import type { WorkspaceImportDeps, NewCompoundJob, NewGeneClusterJob } from "./types";
import type { Session, SessionItem, CompoundItem, GeneClusterItem } from "../session/types";
import { saveSession } from "../session/api";
import { runQuery } from "../query/api";
import { z } from "zod";

export const MAX_ITEMS = 200;

const SubmitJobRespSchema = z.object({
  ok: z.boolean(),
  elapsed_ms: z.number().int().nonnegative(),
  status: z.string().optional(),
}).partial() // we don't actually use the response body here

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
    },
    SubmitJobRespSchema
  )
}

export async function submitGeneClusterJob(
  sessionId: string,
  item: GeneClusterItem,
): Promise<void> {
  await postJson(
    "/api/submitGeneCluster",
    {
      sessionId,
      itemId: item.id,
      name: item.name,
      fileContent: item.fileContent,
    },
    SubmitJobRespSchema
  )
}

// General helper to add items to the session with capacity checks
async function addItemsToSession(
  deps: WorkspaceImportDeps,
  buildItems: (prev: Session, remainingSlots: number) => { updated: Session; newItems: SessionItem[] }
): Promise<{ nextSession: Session | null; newItems: SessionItem[] }> {
  // Load the current session
  const { setSession, pushNotification } = deps;

  let nextSession: Session | null = null;
  let newItems: SessionItem[] = [];

  setSession(prev => {
    const existingCount = prev.items.length;
    const remainingSlots = MAX_ITEMS - existingCount;

    if (remainingSlots <= 0) {
      pushNotification(`Workspace is full. Maximum of ${MAX_ITEMS} items reached.`, "warning");
      nextSession = prev;
      newItems = [];
      return prev;
    }

    const { updated, newItems: createdItems } = buildItems(prev, remainingSlots);

    if (createdItems.length === 0) {
      nextSession = prev;
      newItems = [];
      return prev;
    }

    nextSession = updated;
    newItems = createdItems;
    return updated;
  });

  return { nextSession, newItems };
}

// Submit compounds
export async function importCompoundsBatch(
  deps: WorkspaceImportDeps,
  compounds: NewCompoundJob[]
): Promise<SessionItem[]> {
  const { pushNotification, setSession, sessionId } = deps;

  if (!compounds.length) {
    pushNotification("No valid compounds to import", "warning");
    return [];
  }

  const { nextSession, newItems } = await addItemsToSession(deps, (prev, remainingSlots) => {
    const limitedCompounds =
      compounds.length > remainingSlots ? compounds.slice(0, remainingSlots) : compounds;

    if (limitedCompounds.length < compounds.length) {
      pushNotification(
        `Only ${remainingSlots} compounds were imported due to workspace limit of ${MAX_ITEMS} items.`,
        "warning"
      );
    }

    const createdItems: SessionItem[] = limitedCompounds.map(({ name, smiles }) => ({
      id: crypto.randomUUID(),
      kind: "compound",
      name,
      smiles,
      taggedSmiles: null,
      retrofingerprints: [],
      primarySequences: [],
      status: "queued",
      errorMessage: null,
      updatedAt: Date.now(),
    }));

    return {
      updated: {
        ...prev,
        items: [...prev.items, ...createdItems],
      },
      newItems: createdItems,
    };
  });

  if (!nextSession || newItems.length === 0) {
    // Nothing to do
    return [];
  }

  // Save session before submitting jobs
  try {
    await saveSession(nextSession);
  } catch (err) {
    const msg = err instanceof Error ? err.message : String(err);
    pushNotification(`Failed to save session: ${msg}`, "error");

    const newIds = new Set(newItems.map(ni => ni.id));

    setSession(prev => ({
      ...prev,
      items: prev.items.map(it =>
        newIds.has(it.id)
          ? {
              ...it,
              status: "error",
              errorMessage: `Failed to save session: ${msg}`,
              updatedAt: Date.now(),
            }
          : it
      ),
    }));

    return [];
  }

  // Submit jobs for each new compound (sequential)
  for (const item of newItems) {
    try {
      await submitCompoundJob(sessionId, item as CompoundItem);
    } catch (err) {
      const msg = err instanceof Error ? err.message : String(err);
      pushNotification(`Failed to submit job for compound "${item.name}": ${msg}`, "error");

      setSession(prev => ({
        ...prev,
        items: prev.items.map(it =>
          it.id === item.id
            ? {
                ...it,
                status: "error",
                errorMessage: `Failed to submit job: ${msg}`,
                updatedAt: Date.now(),
              }
            : it
        ),
      }));
    }
  }

  return newItems;
}

// Single compound import wrapper
export async function importCompound(
  deps: WorkspaceImportDeps,
  payload: NewCompoundJob,
): Promise<SessionItem | null> {
  const items = await importCompoundsBatch(deps, [payload]);
  return items[0] ?? null;
}

// Wrapper around importCompound; construct payload by first retrieving compound job info from server
export async function importCompoundById(
  deps: WorkspaceImportDeps,
  compoundId: number,
): Promise<SessionItem | null> {
  const { pushNotification } = deps;

  // Retrieve compound info from server
  let compoundInfo: { name: string; smiles: string };
  try {
    compoundInfo = await runQuery({
      name: "compound_info_by_id",
      params: { "compound_id": compoundId },
      paramSchema: z.object({
        compound_id: z.number().int().positive(),
      }),
    }).then(res => {
      const rows = (res.rows || []) as { name: string; smiles: string }[];
      if (rows.length === 0) {
        throw new Error("No compound found with the given ID");
      }
      return { name: rows[0].name, smiles: rows[0].smiles };
    });
  } catch (err) {
    const msg = err instanceof Error ? err.message : String(err);
    pushNotification(`Failed to retrieve compound info: ${msg}`, "error");
    return null;
  }

  // Import compound
  const item = await importCompound(deps, {
    name: compoundInfo.name,
    smiles: compoundInfo.smiles,
  });

  return item;
}

// Submit gene clusters
export async function importGeneClustersBatch(
  deps: WorkspaceImportDeps,
  clusters: NewGeneClusterJob[]
): Promise<SessionItem[]> {
  const { pushNotification, setSession, sessionId } = deps;

  if (!clusters.length) {
    pushNotification("No gene clusters to import", "warning");
    return [];
  }

  const { nextSession, newItems } = await addItemsToSession(deps, (prev, remainingSlots) => {
    const limited =
      clusters.length > remainingSlots ? clusters.slice(0, remainingSlots) : clusters;

    if (limited.length < clusters.length) {
      pushNotification(
        `Only ${remainingSlots} gene clusters were imported due to workspace limit of ${MAX_ITEMS} items.`,
        "warning"
      );
    }

    const createdItems: SessionItem[] = limited.map(({ name, fileContent }) => ({
      id: crypto.randomUUID(),
      kind: "gene_cluster",
      name,
      fileContent,
      retrofingerprints: [],
      primarySequences: [],
      status: "queued",
      errorMessage: null,
      updatedAt: Date.now(),
    }));

    return {
      updated: {
        ...prev,
        items: [...prev.items, ...createdItems],
      },
      newItems: createdItems,
    };
  });

  if (!nextSession || newItems.length === 0) {
    return [];
  }

  // Save session before submitting jobs
  try {
    await saveSession(nextSession);
  } catch (err) {
    const msg = err instanceof Error ? err.message : String(err);
    pushNotification(`Failed to save session: ${msg}`, "error");

    const newIds = new Set(newItems.map(ni => ni.id));

    setSession(prev => ({
      ...prev,
      items: prev.items.map(it =>
        newIds.has(it.id)
          ? {
              ...it,
              status: "error",
              errorMessage: `Failed to save session: ${msg}`,
              updatedAt: Date.now(),
            }
          : it
      ),
    }));

    return [];
  }

  // Submit jobs for each new gene cluster (sequential)
  for (const item of newItems) {
    try {
      await submitGeneClusterJob(sessionId, item as GeneClusterItem);
    } catch (err) {
      const msg = err instanceof Error ? err.message : String(err);
      pushNotification(`Failed to submit job for gene cluster "${item.name}": ${msg}`, "error");

      setSession(prev => ({
        ...prev,
        items: prev.items.map(it =>
          it.id === item.id
            ? {
                ...it,
                status: "error",
                errorMessage: `Failed to submit job: ${msg}`,
                updatedAt: Date.now(),
              }
            : it
        ),
      }));
    }
  }

  return newItems;
}

// Single gene cluster import wrapper
export async function importGeneCluster(
  deps: WorkspaceImportDeps,
  payload: NewGeneClusterJob,
): Promise<SessionItem | null> {
  const items = await importGeneClustersBatch(deps, [payload]);
  return items[0] ?? null;
}