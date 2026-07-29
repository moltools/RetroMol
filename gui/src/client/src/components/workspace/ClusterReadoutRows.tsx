import React from "react";
import Box from "@mui/material/Box";
import Stack from "@mui/material/Stack";
import Tooltip from "@mui/material/Tooltip";
import Typography from "@mui/material/Typography";
import { horizontalScrollSx } from "../../theme/scrollbarSx";
import type { ClusterModule, ClusterReadout } from "../../features/clusters/types";

// Mirrors the PKS extender-unit classification in retromol_antismash/modules.py
// (PKSModule.substrate): that logic isn't serialized server-side, so it's
// recomputed here from the same anatomy flags for display purposes only.
function pksExtenderLabel(anatomy: Extract<ClusterModule, { type: "PKS" }>["anatomy"]): string {
  const KR = anatomy.has_active_KR;
  const DH = anatomy.has_active_DH;
  const ER = anatomy.has_active_ER;
  const effDH = DH && KR;
  const effER = ER && KR && DH;

  if (!KR) return "PKS_A";
  if (!effDH) return "PKS_B";
  if (!effER) return "PKS_C";
  return "PKS_D";
}

function ModuleChip({ module }: { module: ClusterModule }) {
  const isNRPS = module.type === "NRPS";
  const substrateName = isNRPS ? module.predicted_substrate?.name ?? null : null;
  const hasConfidentSubstrate = isNRPS && !!substrateName && substrateName !== "unknown";

  const label = isNRPS
    ? (hasConfidentSubstrate ? substrateName! : "unknown substrate")
    : pksExtenderLabel(module.anatomy);

  const tooltipLines = [
    `${module.type} module ${module.module_index_in_gene + 1} in ${module.gene_id} (${module.gene_strand} strand)`,
    `domains: ${module.present_domains.join(", ") || "none"}`,
    isNRPS && module.predicted_substrate?.score != null
      ? `prediction confidence: ${(module.predicted_substrate.score * 100).toFixed(1)}%`
      : null,
  ].filter(Boolean).join("\n");

  return (
    <Tooltip title={<span style={{ whiteSpace: "pre-line" }}>{tooltipLines}</span>} arrow>
      <Box
        sx={{
          px: 1.25,
          py: 0.75,
          borderRadius: 1,
          border: "1px solid",
          borderColor: isNRPS
            ? (hasConfidentSubstrate ? "success.main" : "divider")
            : "info.main",
          bgcolor: isNRPS
            ? (hasConfidentSubstrate ? "success.main" : "background.paper")
            : "info.main",
          color: isNRPS
            ? (hasConfidentSubstrate ? "success.contrastText" : "text.primary")
            : "info.contrastText",
          opacity: isNRPS && !hasConfidentSubstrate ? 1 : 0.92,
          fontSize: "0.8rem",
          whiteSpace: "nowrap",
          flexShrink: 0,
          display: "flex",
          flexDirection: "column",
          alignItems: "center",
          lineHeight: 1.3,
          minWidth: 72,
          cursor: "default",
        }}
      >
        <Typography variant="caption" sx={{ opacity: 0.75, fontSize: "0.65rem" }}>
          {module.type} · {module.gene_id}
        </Typography>
        <Typography variant="body2" sx={{ fontWeight: 600 }}>
          {label}
        </Typography>
      </Box>
    </Tooltip>
  );
}

// Read-only view of a gene cluster's parsed linear module readout(s) -- one row
// per antiSMASH region, each a horizontal chain of NRPS/PKS module chips. Unlike
// compounds, there's no "reconstruction" step here: the payload already IS the
// final module chain, so it's rendered directly rather than fetched separately.
export function ClusterReadoutRows({ readouts }: { readouts: ClusterReadout[] }) {
  if (readouts.length === 0) {
    return (
      <Typography variant="body2" color="text.secondary">
        No antiSMASH regions were found in this file.
      </Typography>
    );
  }

  return (
    <Stack spacing={2}>
      {readouts.map((readout, idx) => (
        <Box key={`${readout.id}-${idx}`}>
          <Typography variant="caption" sx={{ color: "text.secondary", fontWeight: 600 }}>
            {readout.id} ({readout.modules.length} module{readout.modules.length === 1 ? "" : "s"})
          </Typography>

          {readout.modules.length === 0 ? (
            <Typography variant="body2" color="text.secondary" sx={{ mt: 0.5 }}>
              No NRPS/PKS modules were detected in this region.
            </Typography>
          ) : (
            <Box sx={{ display: "flex", flexWrap: "nowrap", gap: 1, alignItems: "center", mt: 0.5, ...horizontalScrollSx }}>
              {readout.modules.map((module, mIdx) => (
                <ModuleChip key={`${readout.id}-module-${mIdx}`} module={module} />
              ))}
            </Box>
          )}
        </Box>
      ))}
    </Stack>
  );
}
