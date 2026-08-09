// Shared pairwise/MSA alignment grid rendering, extracted out of WorkspaceDiscovery so
// DialogViewDiscoveryQuery (which renders a *saved* query's results, live or uploaded)
// can reuse the exact same rendering without duplicating it.

import React from "react";
import Box from "@mui/material/Box";
import Button from "@mui/material/Button";
import Chip from "@mui/material/Chip";
import Checkbox from "@mui/material/Checkbox";
import Collapse from "@mui/material/Collapse";
import IconButton from "@mui/material/IconButton";
import Stack from "@mui/material/Stack";
import Typography from "@mui/material/Typography";
import MuiLink from "@mui/material/Link";
import ExpandMoreIcon from "@mui/icons-material/ExpandMore";
import ExpandLessIcon from "@mui/icons-material/ExpandLess";
import DownloadIcon from "@mui/icons-material/Download";
import { useTheme, alpha, type Theme } from "@mui/material/styles";
import type { DiscoveryResult } from "../../features/discovery/types";
import { MotifHoverCard } from "../MotifHoverCard";
import { MotifName } from "../MotifName";
import { horizontalScrollSx } from "../../theme/scrollbarSx";
import { buildAlignmentSvg, downloadSvg, type AlignmentSvgRow } from "./alignmentSvgExport";

// Shared vertical sizing (padding, border width, font size, line height) so every
// row -- label, sequence cells, and score -- resolves to the exact same height.
// If these drift apart, the three columns render as independent stacks whose
// per-row heights don't match, and the mismatch compounds row after row until
// labels/scores are visibly offset from the row they belong to.
export const ROW_CELL_BASE_SX = {
  px: 0.75,
  py: 0.4,
  borderRadius: 1,
  border: "1px dashed transparent",
  fontSize: "0.75rem",
  lineHeight: 1.4,
  boxSizing: "border-box" as const,
};

// Shading intensity range for per-unit similarity highlighting -- kept well short of
// fully transparent/opaque so an occupied cell is always distinguishable from a gap
// (similarity 0) and never drowns out the motif name (similarity 1).
const MATCH_SHADE_MIN_ALPHA = 0.12;
const MATCH_SHADE_MAX_ALPHA = 0.65;

// similarity is a raw Tanimoto value in [0, 1] (identical structure = 1), from the same
// scoring matrix used to compute the alignment -- not rescaled, since 0-1 is already a
// meaningful, comparable scale across every column and every row.
export function similarityColor(theme: Theme, similarity: number | null | undefined): string | undefined {
  if (similarity === null || similarity === undefined) return undefined;
  const clamped = Math.max(0, Math.min(1, similarity));
  return alpha(theme.palette.primary.main, MATCH_SHADE_MIN_ALPHA + clamped * (MATCH_SHADE_MAX_ALPHA - MATCH_SHADE_MIN_ALPHA));
}

function alignedCellSx(name: string | null, columnWidthCh: number, matchColor: string | undefined) {
  return {
    ...ROW_CELL_BASE_SX,
    borderColor: name === null ? "divider" : "primary.main",
    bgcolor: name === null ? "transparent" : matchColor ?? "action.hover",
    width: `${columnWidthCh}ch`,
    textAlign: "center" as const,
    flexShrink: 0,
  };
}

export type AlignmentGridRow = {
  id: string;
  label: string;
  score?: string; // pre-formatted, e.g. "82.3%" -- omit for rows with no score (e.g. the query)
  // Per-column Tanimoto similarity (0-1) against the query's own unit at that same column,
  // null wherever either side is a gap. Index-aligned with `sequence`. Omit entirely for
  // rows with nothing to compare against (e.g. the query row itself), which render with
  // a neutral fill throughout.
  matchStrengths?: (number | null)[];
  sequence: (string | null)[];
};

export function toSvgRows(rows: AlignmentGridRow[]): AlignmentSvgRow[] {
  return rows.map((row) => ({ label: row.label, score: row.score, matchStrengths: row.matchStrengths, sequence: row.sequence }));
}

export function sanitizeFilenamePart(value: string): string {
  return value.replace(/[^a-z0-9]+/gi, "_").replace(/^_+|_+$/g, "") || "alignment";
}

// Renders any number of aligned rows together (2 for a pairwise view, N for an MSA)
// so every column can be sized to fit the longest name across ALL rows at that
// position -- otherwise a short name like "X" lining up under a long one like
// "propionic acid" makes the alignment unreadable. All rows share one scroll
// container so they always scroll in sync, keeping columns aligned.
export function AlignmentGrid({ rows }: { rows: AlignmentGridRow[] }) {
  const theme = useTheme();
  const columnCount = Math.max(0, ...rows.map((r) => r.sequence.length));
  const columnWidths = Array.from({ length: columnCount }, (_, idx) => {
    let width = 1;
    for (const row of rows) {
      width = Math.max(width, row.sequence[idx]?.length ?? 1);
    }
    return width + 1.5; // small buffer beyond the widest label
  });
  const hasScores = rows.some((r) => r.score !== undefined);

  return (
    <Stack direction="row" spacing={0.5}>
      <Stack spacing={0.5} sx={{ flexShrink: 0 }}>
        {rows.map((row) => (
          <Box
            key={row.id}
            sx={{
              ...ROW_CELL_BASE_SX,
              minWidth: 70,
              maxWidth: 160,
              color: "text.secondary",
              fontWeight: 600,
              whiteSpace: "nowrap",
              overflow: "hidden",
              textOverflow: "ellipsis",
            }}
          >
            {row.label}
          </Box>
        ))}
      </Stack>
      <Box sx={{ minWidth: 0, flex: "1 1 0%", ...horizontalScrollSx }}>
        <Stack spacing={0.5}>
          {rows.map((row) => (
            <Stack key={row.id} direction="row" spacing={0.5} sx={{ flexWrap: "nowrap" }}>
              {columnWidths.map((width, idx) => {
                const name = row.sequence[idx] ?? null;
                const matchColor = similarityColor(theme, row.matchStrengths?.[idx]);
                const cell = (
                  <Box sx={alignedCellSx(name, width, matchColor)}>
                    {name === null ? "–" : <MotifName name={name} />}
                  </Box>
                );
                return name === null ? (
                  <React.Fragment key={idx}>{cell}</React.Fragment>
                ) : (
                  <MotifHoverCard key={idx} name={name}>
                    {cell}
                  </MotifHoverCard>
                );
              })}
            </Stack>
          ))}
        </Stack>
      </Box>
      {hasScores && (
        <Stack spacing={0.5} sx={{ flexShrink: 0 }}>
          {rows.map((row) => (
            <Box key={row.id} sx={{ ...ROW_CELL_BASE_SX, minWidth: 56, textAlign: "right" }}>
              {/* a non-breaking space (not "") keeps a text node present so this box's
                  line-height strut matches its siblings' -- an empty string renders no
                  inline content, so the browser omits the strut and the box collapses
                  shorter, throwing off every row below it in this column. */}
              {row.score ?? " "}
            </Box>
          ))}
        </Stack>
      )}
    </Stack>
  );
}

export function DownloadSvgButton({ rows, filename }: { rows: AlignmentGridRow[]; filename: string }) {
  return (
    <Button
      size="small"
      variant="text"
      startIcon={<DownloadIcon fontSize="small" />}
      onClick={() => downloadSvg(buildAlignmentSvg(toSvgRows(rows)), filename)}
    >
      Download SVG
    </Button>
  );
}

export function ResultRow({
  result,
  rank,
  selectedForMsa,
  onToggleSelectedForMsa,
}: {
  result: DiscoveryResult;
  rank: number;
  selectedForMsa: boolean;
  onToggleSelectedForMsa: () => void;
}) {
  const [expanded, setExpanded] = React.useState(false);
  const theme = useTheme();

  return (
    <Box sx={{ borderBottom: "1px solid", borderColor: "divider", py: 1 }}>
      <Stack direction="row" alignItems="center" spacing={1.5}>
        <Checkbox
          size="small"
          checked={selectedForMsa}
          onChange={onToggleSelectedForMsa}
          title="Include in multiple sequence alignment and compound comparison"
          sx={{ p: 0.5 }}
        />

        <Typography variant="body2" sx={{ width: 28, color: "text.secondary" }}>
          {rank}
        </Typography>

        <Box sx={{ flex: 1, minWidth: 0 }}>
          {result.url ? (
            <MuiLink
              href={result.url}
              target="_blank"
              rel="noopener noreferrer"
              underline="hover"
              color={(theme.vars || theme).palette.primary.main}
              sx={{ fontWeight: 500 }}
            >
              {result.name}
            </MuiLink>
          ) : (
            <Typography variant="body2" fontWeight={500} noWrap>
              {result.name}
            </Typography>
          )}
        </Box>

        <Chip
          label={result.origin === "upload" ? "User upload" : result.type === "compound" ? "Compound" : "BGC"}
          size="small"
          color={result.origin === "upload" ? "info" : "default"}
          sx={{ fontSize: "0.7rem" }}
        />

        <Typography variant="caption" sx={{ minWidth: 120, textAlign: "right" }}>
          fingerprint {(result.fingerprintSimilarity * 100).toFixed(1)}%
        </Typography>

        <Typography variant="caption" sx={{ minWidth: 90, textAlign: "right" }}>
          align {result.normalizedAlignmentScorePct.toFixed(1)}%
        </Typography>

        <IconButton size="small" onClick={() => setExpanded((e) => !e)}>
          {expanded ? <ExpandLessIcon fontSize="small" /> : <ExpandMoreIcon fontSize="small" />}
        </IconButton>
      </Stack>

      <Collapse in={expanded} unmountOnExit>
        <Box sx={{ mt: 1, pl: 4.5 }}>
          <AlignmentGrid
            rows={[
              { id: "query", label: "Query", sequence: result.alignedQuery },
              {
                id: result.entryId,
                label: result.type === "compound" ? "Compound" : "BGC",
                score: `${result.normalizedAlignmentScorePct.toFixed(1)}%`,
                matchStrengths: result.alignedSimilarity,
                sequence: result.alignedTarget,
              },
            ]}
          />
          <Box sx={{ mt: 0.5 }}>
            <DownloadSvgButton
              rows={[
                { id: "query", label: "Query", sequence: result.alignedQuery },
                {
                  id: result.entryId,
                  label: result.name,
                  score: `${result.normalizedAlignmentScorePct.toFixed(1)}%`,
                  matchStrengths: result.alignedSimilarity,
                  sequence: result.alignedTarget,
                },
              ]}
              filename={`pairwise_alignment_${sanitizeFilenamePart(result.name)}.svg`}
            />
          </Box>
        </Box>
      </Collapse>
    </Box>
  );
}
