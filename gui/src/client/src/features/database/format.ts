// ChEBI role names (both "chebi_biological_role" and "chebi_chemical_role") come from
// the source data lowercase (e.g. "antibacterial agent", "xenobiotic") -- sentence-cased
// here for display only. Never touches the stored/queried label itself (exact-match
// lookups, external links, and backend grouping all still use the raw lowercase form).
const CHEBI_ROLE_RANKS = new Set(["chebi_biological_role", "chebi_chemical_role"]);

export function isChebiRoleRank(rank: string | null | undefined): boolean {
  return !!rank && CHEBI_ROLE_RANKS.has(rank);
}

export function toSentenceCase(label: string): string {
  // Capitalize only the label's own first character -- "antibacterial agent" ->
  // "Antibacterial agent", not "Antibacterial Agent". Not `\b\w`-regex-based: JS
  // regex `\w` is ASCII-only, so it doesn't recognize letters like "ø" as word
  // characters, which would break `\b` (word-boundary) detection right after one.
  return label.length > 0 ? label.charAt(0).toUpperCase() + label.slice(1) : label;
}
