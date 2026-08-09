import { z } from "zod";

export const MatchingRuleSchema = z.object({
  id: z.string(),
  name: z.string(),
  smiles: z.string(),
  // A friendlier SMILES to depict instead of `smiles`, when set -- `smiles` is written
  // for substructure matching and may carry a reactive leaving-group placeholder.
  displaySmiles: z.string().nullable(),
  pseudonyms: z.array(z.string()),
  terminal: z.boolean(),
  stereochemistry: z.boolean(),
  props: z.record(z.string(), z.any()),
});
export type MatchingRule = z.output<typeof MatchingRuleSchema>;

export const ReactionRuleSchema = z.object({
  id: z.string(),
  name: z.string(),
  smarts: z.string(),
  allowedInBulk: z.boolean(),
  props: z.record(z.string(), z.any()),
});
export type ReactionRule = z.output<typeof ReactionRuleSchema>;

export const RuleSetRespSchema = z.object({
  matchStereochemistry: z.boolean(),
  matchingRules: z.array(MatchingRuleSchema),
  reactionRules: z.array(ReactionRuleSchema),
});
export type RuleSetResp = z.output<typeof RuleSetRespSchema>;
