import { postJson } from "../http";
import { PrimarySequence } from "../session/types";
import { ItemDrawingResultSchema } from "./types";

export async function drawCompoundItem(
  taggedParentSmiles: string,
  primarySequence: PrimarySequence
): Promise<string> {
  const data = await postJson(
    "/api/drawCompoundItem",
    {
      taggedParentSmiles,
      primarySequence
    },
    ItemDrawingResultSchema
  )
  return data.svg;
}

export async function drawGeneClusterItem(
  fileContent: string
): Promise<string> {
  const data = await postJson(
    "/api/drawGeneClusterItem",
    {
      fileContent
    },
    ItemDrawingResultSchema
  )
  return data.svg;
}