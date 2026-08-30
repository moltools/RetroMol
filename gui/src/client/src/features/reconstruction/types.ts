import { z } from "zod";

// Joins two candidate primary sequence paths that couldn't be threaded into a
// single path (e.g. a disconnected sugar, or a branched/cyclic assembly) so a
// compound/BGC always has exactly one primary sequence -- see
// retromol.model.readout.merge_named_paths on the backend. Alignable like any
// other token (see gui/src/server/routes/discovery.py's DiscoveryContext), so it
// can appear in an aligned sequence too, not just an unaligned primary sequence.
export const LINK_TOKEN = "<LINK>";
export const isLinkToken = (name: string | null | undefined): boolean => name === LINK_TOKEN;

// Splits a merged primary sequence back apart on LINK_TOKEN, e.g. so the discovery
// query editor can offer "search with just this chain" alongside "search with the
// whole thing" for a compound whose readout has more than one path (a branched/
// disconnected assembly, or a main chain plus tailoring events like glycosylation).
// Names are [name, tags] pairs (see PrimarySequenceItemSchema) or plain strings --
// either way the link token itself is dropped, never included in a subsequence.
export function splitOnLinkToken<T extends string | PrimarySequenceItem>(sequence: T[]): T[][] {
  const nameOf = (item: T): string => (Array.isArray(item) ? item[0] : item);

  const groups: T[][] = [[]];
  for (const item of sequence) {
    if (isLinkToken(nameOf(item))) {
      groups.push([]);
    } else {
      groups[groups.length - 1].push(item);
    }
  }

  return groups.filter((g) => g.length > 0);
}

export const PrimarySequenceItemSchema = z.tuple([z.string(), z.array(z.number())]);
export type PrimarySequenceItem = z.output<typeof PrimarySequenceItemSchema>;

export const ReconstructionSchema = z.object({
  tagged_input_smiles: z.string(),
  // Null when RetroMol's hardcoded fusion chemistry couldn't rebuild a backbone
  // for this path -- see `backbone_warning` in that case. The primary sequence is
  // still meaningful (it's derived directly from the source structure), so a
  // candidate is kept even without a backbone.
  tagged_backbone_smiles: z.string().nullable(),
  primary_sequence: z.array(PrimarySequenceItemSchema),
  backbone_warning: z.string().nullable().default(null),
  // False when `primary_sequence` isn't a genuine biosynthetic order -- RetroMol
  // identified every building block individually but couldn't connect them into a
  // single path (e.g. a branched or cyclic assembly). Render as an unordered set,
  // not a sequence.
  ordered: z.boolean().default(true),
});
export type Reconstruction = z.output<typeof ReconstructionSchema>;

export const GeneratedBackboneSchema = ReconstructionSchema.extend({
  // Same structure as tagged_backbone_smiles but with RetroMol's internal atom
  // tags stripped -- what a "copy SMILES" affordance should hand the user,
  // rather than leaking isotope tags into a SMILES they'd paste elsewhere. Null
  // under the same conditions tagged_backbone_smiles is.
  backboneSmiles: z.string().nullable().default(null),
});
export type GeneratedBackbone = z.output<typeof GeneratedBackboneSchema>;

export const GenerateBackboneRespSchema = z.object({
  data: GeneratedBackboneSchema,
});

export const ReconstructCompoundRespSchema = z.object({
  ok: z.boolean().optional(),
  status: z.string().optional(),
  data: z.array(ReconstructionSchema).default([]),
  // The single primary sequence read directly off result.linear_readout -- the
  // same representation the database is actually populated with, unlike `data`
  // above (reconstruct_linear_readout's candidates, which apply backbone-
  // reconstruction eligibility/orientation filtering that can diverge from what's
  // stored). Always at most one element. Nothing is filtered out of it -- every
  // path found in the assembly graph, tailoring events (glycosylation,
  // methylation, ...) included, is merged into it, joined by LINK_TOKEN -- so use
  // splitOnLinkToken to pull out just one biosynthetic chain if that's what's
  // wanted for a query, rather than the whole molecule. Query with this, not
  // `data`, to get a fingerprint comparable to the database.
  dbMatchingSequences: z.array(ReconstructionSchema).default([]),
});
