import { postJson } from "../http";
import { ReconstructCompoundRespSchema, Reconstruction } from "./types";

export async function reconstructCompound(
  sessionId: string,
  itemId: string,
  signal?: AbortSignal
): Promise<Reconstruction[]> {
  const data = await postJson(
    "/api/reconstructCompound",
    { sessionId, itemId },
    ReconstructCompoundRespSchema,
    signal
  );
  return data.data;
}
