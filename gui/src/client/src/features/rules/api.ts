import { getJson } from "../http";
import { RuleSetRespSchema, type RuleSetResp } from "./types";
import { ItemDrawingResultSchema } from "../drawing/types";

// The whole default rule set, fetched once and cached by the caller (see
// WorkspaceRules) -- it's small (a few hundred rules total) and effectively static
// for the life of the server process, same rationale as fetchMotifStructures.
export async function fetchRuleSet(signal?: AbortSignal): Promise<RuleSetResp> {
  return getJson("/api/ruleSet", RuleSetRespSchema, signal);
}

// A reaction rule's scheme (reactant(s) -> product(s)), rendered server-side by
// RDKit -- see ReactionScheme for why this replaced client-side rendering.
export async function fetchReactionSchemeSvg(
  ruleId: string,
  theme: "light" | "dark",
  signal?: AbortSignal
): Promise<string> {
  const data = await getJson(
    `/api/reactionSchemeSvg/${encodeURIComponent(ruleId)}?theme=${theme}`,
    ItemDrawingResultSchema,
    signal
  );
  return data.svg;
}
