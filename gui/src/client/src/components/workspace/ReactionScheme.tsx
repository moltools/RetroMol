// Renders a reaction rule's scheme -- reactant(s) -> product(s) -- as an SVG
// rendered server-side by RDKit (see fetchReactionSchemeSvg / routes/rules.py's
// reaction_scheme_svg). This replaced an earlier client-side attempt built on
// smiles-drawer's ReactionDrawer/parseReaction: that API parses with the same
// strict, SMILES-only grammar `SmilesDrawer.parse` uses, so a single SMARTS-only
// token anywhere in the reaction (an OR-list like "[C,c]", a boolean primitive like
// "[O&D2]"/"[C;!R]") aborted the whole parse -- most of this ruleset's reaction
// rules use at least one such token. RDKit's own reaction drawer understands full
// SMARTS query syntax directly, so every rule renders, not just the plain-enough
// handful a SMILES grammar could parse.
import React from "react";
import Alert from "@mui/material/Alert";
import Box from "@mui/material/Box";
import CircularProgress from "@mui/material/CircularProgress";
import DOMPurify from "dompurify";
import { useColorScheme } from "@mui/material/styles";
import { useQuery } from "@tanstack/react-query";
import { fetchReactionSchemeSvg } from "../../features/rules/api";
import { DrawingAttribution } from "../DrawingAttribution";

export function ReactionScheme({ id, smarts }: { id: string; smarts: string }) {
  const { mode, systemMode } = useColorScheme();
  const theme = (systemMode !== undefined ? systemMode : mode) === "dark" ? "dark" : "light";

  const svgQuery = useQuery({
    queryKey: ["reactionSchemeSvg", id, theme],
    queryFn: ({ signal }) => fetchReactionSchemeSvg(id, theme, signal),
    staleTime: Infinity,
    gcTime: Infinity,
  });

  // Server-rendered SVG markup gets injected raw via dangerouslySetInnerHTML, so it
  // must be sanitized first -- same rationale/profile as SvgViewer.
  const sanitizedSvg = React.useMemo(
    () => (svgQuery.data ? DOMPurify.sanitize(svgQuery.data.svg, { USE_PROFILES: { svg: true, svgFilters: true } }) : null),
    [svgQuery.data]
  );

  if (svgQuery.isLoading) {
    return (
      <Box sx={{ display: "flex", alignItems: "center", justifyContent: "center", height: 100 }}>
        <CircularProgress size={20} />
      </Box>
    );
  }

  if (svgQuery.error || !sanitizedSvg || !svgQuery.data) {
    return (
      <Alert severity="warning" variant="outlined" sx={{ py: 0 }}>
        Could not render this reaction ({smarts}).
      </Alert>
    );
  }

  return (
    <>
      <Box
        sx={{ "& svg": { display: "block", maxWidth: "100%", height: "auto" } }}
        dangerouslySetInnerHTML={{ __html: sanitizedSvg }}
      />
      <DrawingAttribution library="rdkit" version={svgQuery.data.rdkitVersion} />
    </>
  );
}
