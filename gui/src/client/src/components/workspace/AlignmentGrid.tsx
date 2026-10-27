// Shared pairwise/MSA alignment grid rendering, extracted out of WorkspaceDiscovery so
// DialogViewDiscoveryQuery (which renders a *saved* query's results, live or uploaded)
// can reuse the exact same rendering without duplicating it.

import React from "react";
import Box from "@mui/material/Box";
import Button from "@mui/material/Button";
import Chip from "@mui/material/Chip";
import Checkbox from "@mui/material/Checkbox";
import CircularProgress from "@mui/material/CircularProgress";
import Collapse from "@mui/material/Collapse";
import Stack from "@mui/material/Stack";
import Tooltip from "@mui/material/Tooltip";
import Typography from "@mui/material/Typography";
import ExpandMoreIcon from "@mui/icons-material/ExpandMore";
import ExpandLessIcon from "@mui/icons-material/ExpandLess";
import DownloadIcon from "@mui/icons-material/Download";
import SendIcon from "@mui/icons-material/Send";
import { useTheme, alpha, type Theme } from "@mui/material/styles";
import type { DiscoveryResult } from "../../features/discovery/types";
import { MotifHoverCard } from "../MotifHoverCard";
import { MotifName } from "../MotifName";
import { horizontalScrollSx } from "../../theme/scrollbarSx";
import { buildAlignmentSvg, downloadSvg, type AlignmentSvgRow } from "./alignmentSvgExport";
import { MinimalIconButton} from "../MinimalIconButton";
import { getEntryAnnotations } from "../../features/database/api";
import type { EntryAnnotation } from "../../features/database/types";
import { isChebiRoleRank, toSentenceCase } from "../../features/database/format";

const ANNOTATION_CATEGORY_LABELS: Record<string, string> = {
  phylogeny: "Phylogeny",
  biosynthetic_class: "Biosynthetic class",
  chemical_class: "Chemical class",
  bioactivity: "Bioactivity",
};

function groupAnnotationsByCategory(annotations: EntryAnnotation[]): [string, EntryAnnotation[]][] {
  const groups = new Map<string, EntryAnnotation[]>();
  for (const a of annotations) {
    const bucket = groups.get(a.category) ?? [];
    bucket.push(a);
    groups.set(a.category, bucket);
  }
  return Array.from(groups.entries());
}

// Shared by the "Databases" (source) row and every AnnotationChips category row --
// one component, not parallel copies of the same Stack/Typography markup, so the
// label-to-chip gap is guaranteed identical by construction rather than by keeping two
// separate JSX blocks in sync by hand.
const AnnotationRow: React.FC<{ label: string; sx?: object; children: React.ReactNode }> = ({
  label,
  sx,
  children,
}) => (
  <Stack direction="row" alignItems="center" sx={{ flexWrap: "wrap", gap: 0.25, ...sx }}>
    <Typography variant="caption" color="text.secondary" sx={{ minWidth: 110 }}>
      {label}
    </Typography>
    {children}
  </Stack>
);

const AnnotationChips: React.FC<{ entryId: string }> = ({ entryId }) => {
  const [annotations, setAnnotations] = React.useState<EntryAnnotation[] | null>(null);
  const [error, setError] = React.useState<string | null>(null);

  React.useEffect(() => {
    const controller = new AbortController();
    getEntryAnnotations(entryId, controller.signal)
      .then((resp) => setAnnotations(resp.results))
      .catch((err) => {
        if (controller.signal.aborted) return;
        setError(err instanceof Error ? err.message : String(err));
      });
    return () => controller.abort();
  }, [entryId]);

  if (error) {
    return (
      <Typography variant="caption" color="text.secondary" sx={{ display: "block", mb: 1 }}>
        Couldn't load annotations: {error}
      </Typography>
    );
  }

  if (annotations === null) {
    return null;
  }

  if (annotations.length === 0) {
    return null;
  }

  return (
    <Stack spacing={0.5} sx={{ mb: 2.5 }}>
      {groupAnnotationsByCategory(annotations).map(([category, items]) => (
        <AnnotationRow key={category} label={ANNOTATION_CATEGORY_LABELS[category] ?? category}>
          {items.map((item) => {
            const label = isChebiRoleRank(item.rank) ? toSentenceCase(item.label) : item.label;
            return item.url ? (
              <Chip
                key={item.id}
                component="a"
                href={item.url}
                target="_blank"
                rel="noopener noreferrer"
                clickable
                color="primary"
                size="small"
                variant="outlined"
                label={label}
              />
            ) : (
              <Chip key={item.id} size="small" variant="outlined" label={label} />
            );
          })}
        </AnnotationRow>
      ))}
    </Stack>
  );
};

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
  const isGap = name === null;
  return {
    ...ROW_CELL_BASE_SX,
    borderColor: isGap ? "divider" : "primary.main",
    // Every cell, gap or occupied, gets an opaque background.paper fill so the
    // row's connecting line (see ROW_LINE_SX) never shows through and breaks
    // the pill's outline. Occupied cells layer the (semi-transparent) match
    // tint as a backgroundImage on top, which composites down to a fully
    // solid color over the opaque base.
    bgcolor: "background.paper",
    backgroundImage: !isGap && matchColor ? `linear-gradient(${matchColor}, ${matchColor})` : "none",
    width: `${columnWidthCh}ch`,
    textAlign: "center" as const,
    flexShrink: 0,
    // Sits above the row's connecting line (see ROW_LINE_SX).
    position: "relative" as const,
    zIndex: 1,
  };
}

// The "string" that runs behind one row of aligned units -- the same trick
// PrimarySequenceChips uses for a sequence of chips -- so the eye can follow a
// row across its gaps without losing the thread. Deliberately NOT extended to
// the label/score columns: those are plain, variable-width text (not a fixed
// row of bordered chips), so hiding the line behind them would need an opaque
// backdrop wide enough to cover arbitrarily long/short labels, which either
// shows as a mismatched-looking box around the text or (worse, if given its
// own wrapper element) risks throwing off the shared row height those columns
// depend on to stay aligned with the sequence column beside them.
const ROW_LINE_SX = {
  position: "relative" as const,
  "&::before": {
    content: '""',
    position: "absolute" as const,
    left: 0,
    right: 0,
    top: "50%",
    height: "2px",
    bgcolor: "divider",
    zIndex: 0,
  },
};

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
            <Stack key={row.id} direction="row" spacing={0.5} sx={{ flexWrap: "nowrap", ...ROW_LINE_SX }}>
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
              {/* A non-breaking space -- written as an escape, not a literal character, since a
                  literal nbsp glyph in source is invisible and indistinguishable from a plain
                  space in an editor, exactly how this regressed silently before: a plain " "
                  collapses under normal whitespace rules (it has nothing non-whitespace to
                  anchor to), so the box loses its line-height strut and shrinks, throwing every
                  row below it in this column out of alignment with the label/sequence columns. */}
              {row.score ?? "\u00A0"}
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
  onSendToUploads,
  uploadSlotsFull,
  sendingToUploads,
}: {
  result: DiscoveryResult;
  rank: number;
  selectedForMsa: boolean;
  onToggleSelectedForMsa: () => void;
  // Only offered for database-origin results -- an "upload" origin result is already
  // sitting in the Uploads tab, so sending it back would just duplicate it. Omit
  // entirely (rather than passing a no-op) when the caller has no session to add to.
  onSendToUploads?: (result: DiscoveryResult) => void;
  // Whether the Uploads tab is already at MAX_ITEMS -- disables the button with an
  // explanatory tooltip instead of letting the user find out only after clicking.
  uploadSlotsFull?: boolean;
  sendingToUploads?: boolean;
}) {
  const [expanded, setExpanded] = React.useState(false);

  const canOfferSend = result.origin === "database" && !!onSendToUploads;
  const hasRaw = !!result.raw;
  const sendDisabledReason = !hasRaw
    ? `No original ${result.type === "compound" ? "SMILES" : "GenBank"} data was stored for this entry, so it can't be reparsed.`
    : uploadSlotsFull
    ? "Uploads tab is already at its maximum of 50 items."
    : null;

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
          <Typography variant="body2" fontWeight={500} noWrap>
            {result.name}
          </Typography>
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

        {canOfferSend && (
          <Tooltip title={sendDisabledReason ?? "Send to Uploads to reparse with the current rule set"} arrow>
            <span>
              <MinimalIconButton
                size="small"
                disabled={!!sendDisabledReason || sendingToUploads}
                onClick={() => onSendToUploads!(result)}
              >
                {sendingToUploads ? <CircularProgress size={16} /> : <SendIcon fontSize="small" />}
              </MinimalIconButton>
            </span>
          </Tooltip>
        )}

        <MinimalIconButton size="small" onClick={() => setExpanded((e) => !e)}>
          {expanded ? <ExpandLessIcon fontSize="small" /> : <ExpandMoreIcon fontSize="small" />}
        </MinimalIconButton>
      </Stack>

      <Collapse in={expanded} unmountOnExit>
        <Box sx={{ mt: 1, pl: 4.5 }}>
          {result.sources.length > 0 && (
            <AnnotationRow label="Databases" sx={{ mb: 0.5 }}>
              {result.sources.map((source, idx) =>
                source.url ? (
                  <Chip
                    key={idx}
                    component="a"
                    href={source.url}
                    target="_blank"
                    rel="noopener noreferrer"
                    clickable
                    color="primary"
                    size="small"
                    variant="outlined"
                    label={`${source.databaseName}: ${source.name}`}
                  />
                ) : (
                  <Chip
                    key={idx}
                    size="small"
                    variant="outlined"
                    label={`${source.databaseName}: ${source.name}`}
                  />
                )
              )}
            </AnnotationRow>
          )}
          <AnnotationChips entryId={result.entryId} />
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
