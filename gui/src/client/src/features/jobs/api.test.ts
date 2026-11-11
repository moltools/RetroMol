import { importCompoundsBatch, MAX_ITEMS } from "./api";
import * as sessionApi from "../session/api";
import * as http from "../http";
import type { Session, SessionItem } from "../session/types";
import type { WorkspaceImportDeps } from "./types";

jest.mock("../session/api");
jest.mock("../http");

const mockSaveSession = sessionApi.saveSession as jest.Mock;
const mockPostJson = http.postJson as jest.Mock;

function makeExistingCompound(id: string): SessionItem {
  return {
    id,
    kind: "compound",
    name: id,
    smiles: "C",
    matchStereochemistry: false,
    status: "done",
    errorMessage: null,
    updatedAt: Date.now(),
    score: 1,
  };
}

// Drives the same (session, setSession) contract WorkspaceUpload builds for
// real, but backed by a plain local variable instead of React state.
function makeDeps(initialSession: Session) {
  let session = initialSession;
  const pushNotification = jest.fn();
  const setSession = (updater: (prev: Session) => Session) => {
    session = updater(session);
  };
  const deps: WorkspaceImportDeps = { setSession, pushNotification, sessionId: session.sessionId };
  return { deps, pushNotification, getSession: () => session };
}

describe("importCompoundsBatch", () => {
  beforeEach(() => {
    jest.resetAllMocks();
    mockSaveSession.mockResolvedValue(undefined);
    mockPostJson.mockResolvedValue({ ok: true });
  });

  it("caps the batch at MAX_ITEMS and warns about anything dropped", async () => {
    const existing = Array.from(
      { length: MAX_ITEMS - 2 },
      (_, i) => makeExistingCompound(`existing-${i}`)
    );
    const { deps, getSession, pushNotification } = makeDeps({
      sessionId: "s1",
      created: Date.now(),
      items: existing,
    });

    const compounds = Array.from({ length: 5 }, (_, i) => ({
      name: `new${i}`,
      smiles: "CCO",
      matchStereochemistry: false,
    }));

    const result = await importCompoundsBatch(deps, compounds);

    expect(result).toHaveLength(2); // only 2 slots were left
    expect(getSession().items).toHaveLength(MAX_ITEMS);
    expect(pushNotification).toHaveBeenCalledWith(
      expect.stringContaining("Only importing 2 compounds"),
      "warning"
    );
  });

  it("refuses to import into an already-full session", async () => {
    const existing = Array.from({ length: MAX_ITEMS }, (_, i) => makeExistingCompound(`existing-${i}`));
    const { deps, getSession, pushNotification } = makeDeps({
      sessionId: "s1",
      created: Date.now(),
      items: existing,
    });

    const result = await importCompoundsBatch(deps, [{ name: "n", smiles: "C", matchStereochemistry: false }]);

    expect(result).toHaveLength(0);
    expect(getSession().items).toHaveLength(MAX_ITEMS);
    expect(pushNotification).toHaveBeenCalledWith(
      expect.stringContaining(`maximum of ${MAX_ITEMS} items`),
      "warning"
    );
  });

  it("marks only the failing item as an error, without losing the rest of the batch", async () => {
    mockPostJson.mockImplementation((_url: string, body: any) => {
      if (body?.name === "bad") return Promise.reject(new Error("boom"));
      return Promise.resolve({ ok: true });
    });

    const { deps, getSession, pushNotification } = makeDeps({
      sessionId: "s1",
      created: Date.now(),
      items: [],
    });

    await importCompoundsBatch(deps, [
      { name: "good1", smiles: "C", matchStereochemistry: false },
      { name: "bad", smiles: "C", matchStereochemistry: false },
      { name: "good2", smiles: "C", matchStereochemistry: false },
    ]);

    const items = getSession().items;
    const bad = items.find((it) => it.name === "bad");
    const goods = items.filter((it) => it.name !== "bad");

    expect(bad?.status).toBe("error");
    expect(bad?.errorMessage).toContain("boom");
    goods.forEach((it) => expect(it.status).toBe("queued"));

    expect(pushNotification).toHaveBeenCalledWith(
      expect.stringContaining('Failed to submit job for compound "bad"'),
      "error"
    );
  });

  it("marks the whole batch as errored if the session fails to save", async () => {
    mockSaveSession.mockRejectedValue(new Error("redis down"));

    const { deps, getSession, pushNotification } = makeDeps({
      sessionId: "s1",
      created: Date.now(),
      items: [],
    });

    const result = await importCompoundsBatch(deps, [{ name: "n", smiles: "C", matchStereochemistry: false }]);

    expect(result).toHaveLength(0);
    expect(getSession().items[0].status).toBe("error");
    expect(mockPostJson).not.toHaveBeenCalled(); // never got to submitting jobs
    expect(pushNotification).toHaveBeenCalledWith(
      expect.stringContaining("Failed to save session before importing compounds"),
      "error"
    );
  });
});
