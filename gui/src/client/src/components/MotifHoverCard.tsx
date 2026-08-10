import React from "react";
import Box from "@mui/material/Box";
import CircularProgress from "@mui/material/CircularProgress";
import Tooltip from "@mui/material/Tooltip";
import Typography from "@mui/material/Typography";
import { useQuery } from "@tanstack/react-query";
import { fetchMotifStructures } from "../features/motifs/api";
import { MotifName } from "./MotifName";
import SmilesDrawerContainer from "./SmilesDrawerContainer.js";
import { DrawingAttribution } from "./DrawingAttribution";

const DRAWING_SIZE = 100;

// Content of the popup -- a separate component (rather than inlined in the
// Tooltip's title) so its useQuery only actually runs once the Tooltip mounts it,
// which MUI's Popper does lazily on open. queryKey is shared across every
// instance on the page, so the name -> SMILES map is fetched once and cached.
function MotifHoverContent({ name, hint }: { name: string; hint?: string }) {
  const structuresQuery = useQuery({
    queryKey: ["motifStructures"],
    queryFn: ({ signal }) => fetchMotifStructures(signal),
    staleTime: Infinity,
    gcTime: Infinity,
  });

  const smiles = structuresQuery.data?.[name];

  // Must be unique per *mounted instance*, not derived from `name` alone --
  // SmilesDrawerContainer draws into a DOM node it looks up by this id, and the
  // same motif name commonly appears in many cells at once (a whole column of
  // an MSA, a run of identical residues). Two same-named hover cards can end up
  // mounted at the same moment -- e.g. dragging the pointer straight down a
  // repeated column, where the outgoing tooltip's exit transition briefly
  // overlaps the incoming one's mount -- and a shared id means the lookup can
  // resolve to the wrong (closing) node, leaving the one you're actually
  // hovering empty until it flashes into view as the other one unmounts.
  const reactId = React.useId();

  return (
    <Box sx={{ display: "flex", flexDirection: "column", alignItems: "center", gap: 0.5, maxWidth: 200, py: 0.5 }}>
      {smiles ? (
        <>
          <SmilesDrawerContainer identifier={`motif-hover-${reactId}`} smiles={smiles} size={DRAWING_SIZE} />
          <DrawingAttribution library="smiles-drawer" sx={{ fontSize: "0.65rem" }} />
        </>
      ) : (
        <Box
          sx={{
            width: DRAWING_SIZE,
            height: DRAWING_SIZE,
            display: "flex",
            alignItems: "center",
            justifyContent: "center",
          }}
        >
          {structuresQuery.isLoading ? (
            <CircularProgress size={20} />
          ) : (
            <Typography variant="caption" color="text.secondary" sx={{ textAlign: "center" }}>
              No structure available
            </Typography>
          )}
        </Box>
      )}
      <Typography variant="caption" fontWeight={600}>
        <MotifName name={name} />
      </Typography>
      {hint && (
        <Typography variant="caption" color="text.secondary" sx={{ textAlign: "center" }}>
          {hint}
        </Typography>
      )}
    </Box>
  );
}

// Wraps a motif chip/cell (primary sequence chips, the sequence editor's blocks,
// and the shared pairwise/MSA AlignmentGrid all render a bare motif name) so
// hovering it shows a small structure preview, its name, and optionally whatever
// contextual hint the caller already showed as a plain-text tooltip.
export function MotifHoverCard({
  name,
  hint,
  children,
}: {
  name: string;
  hint?: string;
  children: React.ReactElement;
}) {
  return (
    <Tooltip
      title={<MotifHoverContent name={name} hint={hint} />}
      arrow
      enterDelay={400}
      slotProps={{
        tooltip: {
          sx: {
            bgcolor: "background.paper",
            color: "text.primary",
            border: "1px solid",
            borderColor: "divider",
            boxShadow: 3,
          },
        },
        arrow: {
          sx: {
            color: "background.paper",
            "&::before": {
              border: "1px solid",
              borderColor: "divider",
            },
          },
        },
      }}
    >
      {children}
    </Tooltip>
  );
}
