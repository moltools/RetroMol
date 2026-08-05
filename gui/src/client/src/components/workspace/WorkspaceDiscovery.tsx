import React from "react";
import Box from "@mui/material/Box";
import Card from "@mui/material/Card";
import CardContent from "@mui/material/CardContent";
import Stack from "@mui/material/Stack";
import Typography from "@mui/material/Typography";
import Alert from "@mui/material/Alert";
import Button from "@mui/material/Button";
import Checkbox from "@mui/material/Checkbox";
import Chip from "@mui/material/Chip";
import CircularProgress from "@mui/material/CircularProgress";
import Collapse from "@mui/material/Collapse";
import Divider from "@mui/material/Divider";
import FormControlLabel from "@mui/material/FormControlLabel";
import ListSubheader from "@mui/material/ListSubheader";
import MenuItem from "@mui/material/MenuItem";
import TextField from "@mui/material/TextField";
import ToggleButton from "@mui/material/ToggleButton";
import ToggleButtonGroup from "@mui/material/ToggleButtonGroup";
import MuiLink from "@mui/material/Link";
import IconButton from "@mui/material/IconButton";
import ExpandMoreIcon from "@mui/icons-material/ExpandMore";
import ExpandLessIcon from "@mui/icons-material/ExpandLess";
import DownloadIcon from "@mui/icons-material/Download";
import { useQuery, useMutation } from "@tanstack/react-query";
import { useTheme, alpha, type Theme } from "@mui/material/styles";
import { Session, SessionItem } from "../../features/session/types";
import { reconstructCompound } from "../../features/reconstruction/api";
import type { PrimarySequenceItem } from "../../features/reconstruction/types";
import { reconstructGeneCluster } from "../../features/clusters/api";
import type { ClusterPrimarySequence } from "../../features/clusters/types";
import { runDiscoveryQuery, runDiscoveryMsa } from "../../features/discovery/api";
import {
  DISCOVERY_SCORE_MODE_OPTIONS,
  type DiscoveryEntryType,
  type DiscoveryResult,
  type DiscoveryScoreMode,
} from "../../features/discovery/types";
import { useNotifications } from "../NotificationProvider";
import { MotifName } from "../MotifName";
import { horizontalScrollSx } from "../../theme/scrollbarSx";
import { buildAlignmentSvg, downloadSvg, type AlignmentSvgRow } from "./alignmentSvgExport";
import { SequenceEditor, type SequenceBlock } from "./SequenceEditor";
import { WorkspaceCompare } from "./WorkspaceCompare";
import { WorkspaceShape } from "./WorkspaceShape";

type WorkspaceDiscoveryProps = {
  session: Session;
  setSession: React.Dispatch<React.SetStateAction<Session | null>>;
};

const MAX_TOP_X = 200;

function blocksFromNames(names: string[]): SequenceBlock[] {
  return names.map((name) => ({ id: crypto.randomUUID(), name }));
}

// Read-only preview of a candidate reconstruction path, for picking which one seeds the editor.
function ReconstructionPreview({ sequence }: { sequence: PrimarySequenceItem[] }) {
  return <NamesPreview names={sequence.map(([name]) => name)} />;
}

// Read-only preview of one BGC region's primary sequence, for picking which one seeds
// the editor. Unlike a compound reconstruction there are no per-block tags to carry,
// just names -- see retromol_antismash.modules.bgc_primary_sequence.
function NamesPreview({ names }: { names: string[] }) {
  return (
    <Box sx={{ display: "flex", flexWrap: "nowrap", gap: 0.5, minWidth: 0, flex: "1 1 0%", ...horizontalScrollSx }}>
      {names.map((name, idx) => (
        <Box
          key={`${name}-${idx}`}
          sx={{
            px: 1,
            py: 0.5,
            borderRadius: 1,
            border: "1px solid",
            borderColor: "divider",
            fontSize: "0.8rem",
            whiteSpace: "nowrap",
            flexShrink: 0,
          }}
        >
          <MotifName name={name} />
        </Box>
      ))}
    </Box>
  );
}

// Shared vertical sizing (padding, border width, font size, line height) so every
// row -- label, sequence cells, and score -- resolves to the exact same height.
// If these drift apart, the three columns render as independent stacks whose
// per-row heights don't match, and the mismatch compounds row after row until
// labels/scores are visibly offset from the row they belong to.
const ROW_CELL_BASE_SX = {
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
function similarityColor(theme: Theme, similarity: number | null | undefined): string | undefined {
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
  // Per-column Tanimoto similarity (0-1) against the query's unit at that same column,
  // null wherever either side is a gap. Index-aligned with `sequence`. Omit entirely for
  // rows with nothing to compare against (e.g. the query row itself), which render with
  // a neutral fill throughout.
  matchStrengths?: (number | null)[];
  sequence: (string | null)[];
};

function toSvgRows(rows: AlignmentGridRow[]): AlignmentSvgRow[] {
  return rows.map((row) => ({ label: row.label, score: row.score, matchStrengths: row.matchStrengths, sequence: row.sequence }));
}

function sanitizeFilenamePart(value: string): string {
  return value.replace(/[^a-z0-9]+/gi, "_").replace(/^_+|_+$/g, "") || "alignment";
}

// Renders any number of aligned rows together (2 for a pairwise view, N for an MSA)
// so every column can be sized to fit the longest name across ALL rows at that
// position -- otherwise a short name like "X" lining up under a long one like
// "propionic acid" makes the alignment unreadable. All rows share one scroll
// container so they always scroll in sync, keeping columns aligned.
function AlignmentGrid({ rows }: { rows: AlignmentGridRow[] }) {
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
                return (
                  <Box key={idx} sx={alignedCellSx(name, width, matchColor)}>
                    {name === null ? "–" : <MotifName name={name} />}
                  </Box>
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
              {row.score ?? " "}
            </Box>
          ))}
        </Stack>
      )}
    </Stack>
  );
}

function DownloadSvgButton({ rows, filename }: { rows: AlignmentGridRow[]; filename: string }) {
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

function ResultRow({
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

        <Typography variant="caption" sx={{ minWidth: 90, textAlign: "right" }}>
          fp {(result.fingerprintSimilarity * 100).toFixed(1)}%
        </Typography>

        <Typography variant="caption" sx={{ minWidth: 100, textAlign: "right" }}>
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

export const WorkspaceDiscovery: React.FC<WorkspaceDiscoveryProps> = ({ session }) => {
  const { pushNotification } = useNotifications();

  const compoundItems = session.items.filter((item): item is SessionItem & { kind: "compound" } => item.kind === "compound");
  const clusterItems = session.items.filter((item): item is SessionItem & { kind: "cluster" } => item.kind === "cluster");

  const [selectedItemId, setSelectedItemId] = React.useState<string>("");
  const selectedItem = session.items.find((item) => item.id === selectedItemId);
  const [blocks, setBlocks] = React.useState<SequenceBlock[]>([]);

  const [entryType, setEntryType] = React.useState<DiscoveryEntryType>("compound");
  const [scoreMode, setScoreMode] = React.useState<DiscoveryScoreMode>("subsequence");
  const [n, setN] = React.useState<number>(100);
  const [topX, setTopX] = React.useState<number>(20);
  const [includeUserUploads, setIncludeUserUploads] = React.useState<boolean>(false);
  const [onlyUserUploads, setOnlyUserUploads] = React.useState<boolean>(false);
  const [resultsView, setResultsView] = React.useState<"pairwise" | "msa" | "compare" | "shape">("pairwise");
  const [selectedForMsa, setSelectedForMsa] = React.useState<Set<string>>(new Set());
  // The originating compound's own SMILES, captured at query time -- needed for the
  // Compare view's structure-based Tanimoto metric, which has nothing to do with a
  // sequence of blocks and can't be recovered from `blocks` after the fact. Null
  // whenever the query wasn't seeded from a compound with a known structure (e.g.
  // built from scratch in the sequence editor, or from a gene cluster).
  const [queryOriginSmiles, setQueryOriginSmiles] = React.useState<string | null>(null);

  const reconstructionQuery = useQuery({
    queryKey: ["reconstructCompound", session.sessionId, selectedItemId],
    queryFn: ({ signal }) => reconstructCompound(session.sessionId, selectedItemId, signal),
    enabled: selectedItem?.kind === "compound",
  });

  const clusterReconstructionQuery = useQuery({
    queryKey: ["reconstructGeneCluster", session.sessionId, selectedItemId],
    queryFn: ({ signal }) => reconstructGeneCluster(session.sessionId, selectedItemId, signal),
    enabled: selectedItem?.kind === "cluster",
  });

  const discoveryMutation = useMutation({
    mutationFn: () =>
      runDiscoveryQuery({
        primarySequence: blocks.map((b) => b.name),
        entryType,
        scoreMode,
        n,
        topX,
        includeUserUploads: includeUserUploads || onlyUserUploads,
        onlyUserUploads,
        sessionId: session.sessionId,
      }),
    onError: (err) => {
      const msg = err instanceof Error ? err.message : String(err);
      pushNotification(`Discovery query failed: ${msg}`, "error");
    },
    onSuccess: (data) => {
      // Fresh results default to "everyone's in" for the MSA; the user can narrow from there.
      setSelectedForMsa(new Set(data.results.map((r) => r.entryId)));
    },
  });

  const toggleSelectedForMsa = React.useCallback((entryId: string) => {
    setSelectedForMsa((prev) => {
      const next = new Set(prev);
      if (next.has(entryId)) next.delete(entryId);
      else next.add(entryId);
      return next;
    });
  }, []);

  const selectedResults = React.useMemo(
    () => discoveryMutation.data?.results.filter((r) => selectedForMsa.has(r.entryId)) ?? [],
    [discoveryMutation.data, selectedForMsa]
  );

  // Anchored on the query as the star center, so every result is ordered/oriented
  // relative to it -- matching how the pairwise view is already anchored. Keyed off
  // the selected members so switching view modes back and forth doesn't refetch, but
  // running a new query or changing the selection naturally does. staleTime is what
  // actually enforces that: without it, toggling `enabled` back to true (i.e.
  // re-entering this view) would refetch regardless of the key being unchanged,
  // since a query is stale-on-mount/re-enable by default (staleTime 0).
  const msaQuery = useQuery({
    queryKey: [
      "discoveryMsa",
      discoveryMutation.data?.querySequence,
      selectedResults.map((r) => r.entryId),
    ],
    queryFn: ({ signal }) =>
      runDiscoveryMsa(
        {
          querySequence: discoveryMutation.data!.querySequence,
          sequences: selectedResults.map((r) => ({ id: r.entryId, sequence: r.primarySequence })),
        },
        signal
      ),
    enabled: resultsView === "msa" && selectedResults.length > 0,
    staleTime: Infinity,
  });

  // Scores aren't recomputed for the MSA -- each row's score is the same
  // normalizedAlignmentScorePct already shown in the pairwise view, re-attached by id,
  // so the two views never disagree on a number (and shade the same way, too).
  const msaRows: AlignmentGridRow[] = React.useMemo(() => {
    if (!msaQuery.data || !discoveryMutation.data) return [];
    const resultsById = new Map(discoveryMutation.data.results.map((r) => [r.entryId, r]));
    return msaQuery.data.rows.map((row) => {
      if (row.id === "query") {
        return { id: "query", label: "Query", sequence: row.alignedSequence };
      }
      const result = resultsById.get(row.id);
      return {
        id: row.id,
        label: result?.name ?? row.id,
        score: result ? `${result.normalizedAlignmentScorePct.toFixed(1)}%` : undefined,
        matchStrengths: row.similarityToQuery,
        sequence: row.alignedSequence,
      };
    });
  }, [msaQuery.data, discoveryMutation.data]);

  const handlePickNames = (names: string[]) => {
    setBlocks(blocksFromNames(names));
  };

  const maxTopX = Math.max(1, Math.min(n, MAX_TOP_X));

  const canQuery = blocks.length > 0 && !discoveryMutation.isPending;

  return (
    <Box sx={{ width: "100%", mx: "auto", display: "flex", flexDirection: "column", gap: "16px" }}>
      <Card variant="outlined">
        <CardContent>
          <Typography component="h1" variant="subtitle1">
            Pick a starting sequence
          </Typography>
          <Typography variant="body2" color="text.secondary" sx={{ mb: 1.5 }}>
            Seed the editor below from one of your uploaded compounds or gene clusters, or build a sequence from scratch by adding blocks directly.
          </Typography>

          {compoundItems.length === 0 && clusterItems.length === 0 ? (
            <Typography variant="body2" color="text.secondary">
              Nothing uploaded yet. Import a compound or gene cluster from the Upload tab first.
            </Typography>
          ) : (
            <TextField
              select
              size="small"
              label="Compound or gene cluster"
              value={selectedItemId}
              onChange={(e) => setSelectedItemId(e.target.value)}
              sx={{ minWidth: 260 }}
            >
              {[
                ...(compoundItems.length > 0 ? [<ListSubheader key="compounds-header">Compounds</ListSubheader>] : []),
                ...compoundItems.map((item) => (
                  <MenuItem key={item.id} value={item.id} disabled={item.status !== "done"}>
                    {item.name} {item.status !== "done" ? `(${item.status})` : ""}
                  </MenuItem>
                )),
                ...(clusterItems.length > 0 ? [<ListSubheader key="clusters-header">Gene clusters</ListSubheader>] : []),
                ...clusterItems.map((item) => (
                  <MenuItem key={item.id} value={item.id} disabled={item.status !== "done"}>
                    {item.name} {item.status !== "done" ? `(${item.status})` : ""}
                  </MenuItem>
                )),
              ]}
            </TextField>
          )}

          {selectedItem?.kind === "compound" && (
            <Box sx={{ mt: 2 }}>
              {reconstructionQuery.isLoading && <CircularProgress size={20} />}
              {reconstructionQuery.error && (
                <Alert severity="error">
                  {(reconstructionQuery.error as Error).message || "Failed to load reconstruction."}
                </Alert>
              )}
              {reconstructionQuery.data && reconstructionQuery.data.length === 0 && (
                <Typography variant="body2" color="text.secondary">
                  No reconstructed primary sequences found for this compound.
                </Typography>
              )}
              <Stack spacing={1}>
                {(reconstructionQuery.data ?? []).map((reconstruction, idx) => {
                  // Prefer whatever was saved for this reconstruction in the Upload
                  // tab's viewer over the raw algorithm output, so a correction made
                  // there is what actually gets queried here.
                  const override = selectedItem.editedPrimarySequences?.[String(idx)];
                  const effectiveSequence = override ?? reconstruction.primary_sequence;

                  return (
                    <Box
                      key={idx}
                      sx={{
                        p: 1,
                        borderRadius: 1,
                        border: "1px solid",
                        borderColor: "divider",
                        display: "flex",
                        alignItems: "center",
                        gap: 1.5,
                        flexWrap: "wrap",
                      }}
                    >
                      <ReconstructionPreview sequence={effectiveSequence} />
                      {override && (
                        <Chip label="Edited" size="small" color="info" variant="outlined" sx={{ fontSize: "0.7rem" }} />
                      )}
                      <Button
                        size="small"
                        variant="outlined"
                        onClick={() => handlePickNames(effectiveSequence.map(([name]) => name))}
                      >
                        Use this
                      </Button>
                    </Box>
                  );
                })}
              </Stack>
            </Box>
          )}

          {selectedItem?.kind === "cluster" && (
            <Box sx={{ mt: 2 }}>
              {clusterReconstructionQuery.isLoading && <CircularProgress size={20} />}
              {clusterReconstructionQuery.error && (
                <Alert severity="error">
                  {(clusterReconstructionQuery.error as Error).message || "Failed to load gene cluster readout."}
                </Alert>
              )}
              {clusterReconstructionQuery.data && clusterReconstructionQuery.data.length === 0 && (
                <Typography variant="body2" color="text.secondary">
                  No antiSMASH regions found for this gene cluster.
                </Typography>
              )}
              <Stack spacing={1}>
                {(clusterReconstructionQuery.data ?? []).map((region: ClusterPrimarySequence) => (
                  <Box
                    key={region.id}
                    sx={{
                      p: 1,
                      borderRadius: 1,
                      border: "1px solid",
                      borderColor: "divider",
                      display: "flex",
                      alignItems: "center",
                      gap: 1.5,
                      flexWrap: "wrap",
                    }}
                  >
                    <Typography variant="caption" color="text.secondary" sx={{ flexShrink: 0 }}>
                      {region.id}
                    </Typography>
                    <NamesPreview names={region.primary_sequence} />
                    <Button size="small" variant="outlined" onClick={() => handlePickNames(region.primary_sequence)}>
                      Use this
                    </Button>
                  </Box>
                ))}
              </Stack>
            </Box>
          )}
        </CardContent>
      </Card>

      <Card variant="outlined">
        <CardContent>
          <Typography component="h1" variant="subtitle1">
            Edit sequence
          </Typography>
          <Typography variant="body2" color="text.secondary" sx={{ mb: 1.5 }}>
            Drag blocks to reorder, click × to delete, or click + to add a new block at the end. You can add stereochemistry-aware motifs by using the '^' symbol.
          </Typography>

          {blocks.length === 0 ? (
            <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
              No blocks yet.
            </Typography>
          ) : null}

          <SequenceEditor blocks={blocks} onChange={setBlocks} disabled={discoveryMutation.isPending} />
        </CardContent>
      </Card>

      <Card variant="outlined">
        <CardContent>
          <Typography component="h1" variant="subtitle1">
            Query
          </Typography>

          <Stack direction="row" spacing={2} alignItems="center" flexWrap="wrap" sx={{ mt: 1.5 }}>
            <ToggleButtonGroup
              size="small"
              exclusive
              value={entryType}
              onChange={(_, value) => value && setEntryType(value)}
              disabled={discoveryMutation.isPending}
            >
              <ToggleButton value="compound">Compounds</ToggleButton>
              <ToggleButton value="bgc">BGCs</ToggleButton>
              <ToggleButton value="both">Both</ToggleButton>
            </ToggleButtonGroup>

            <TextField
              select
              size="small"
              label="Score mode"
              value={scoreMode}
              onChange={(e) => setScoreMode(e.target.value as DiscoveryScoreMode)}
              disabled={discoveryMutation.isPending}
              sx={{ width: 220 }}
            >
              {DISCOVERY_SCORE_MODE_OPTIONS.map((option) => (
                <MenuItem key={option.value} value={option.value}>
                  {option.label}
                </MenuItem>
              ))}
            </TextField>

            <TextField
              size="small"
              type="number"
              label="Retrieve closest N"
              value={n}
              onChange={(e) => {
                const value = Math.max(1, Math.min(1000, Number(e.target.value) || 1));
                setN(value);
                setTopX((prev) => Math.min(prev, Math.max(1, Math.min(value, MAX_TOP_X))));
              }}
              slotProps={{ htmlInput: { min: 1, max: 1000 } }}
              sx={{ width: 160 }}
              disabled={discoveryMutation.isPending}
            />

            <TextField
              size="small"
              type="number"
              label="Show top X"
              value={topX}
              onChange={(e) => setTopX(Math.max(1, Math.min(maxTopX, Number(e.target.value) || 1)))}
              slotProps={{ htmlInput: { min: 1, max: maxTopX } }}
              sx={{ width: 140 }}
              disabled={discoveryMutation.isPending}
            />

            <FormControlLabel
              control={
                <Checkbox
                  size="small"
                  checked={includeUserUploads || onlyUserUploads}
                  onChange={(e) => setIncludeUserUploads(e.target.checked)}
                  disabled={discoveryMutation.isPending || compoundItems.length === 0 || onlyUserUploads}
                />
              }
              label="Include my uploads"
              title="Uploaded compounds and gene clusters compete for a spot among the nearest N candidates, then follow the usual top-X ranking."
            />

            <FormControlLabel
              control={
                <Checkbox
                  size="small"
                  checked={onlyUserUploads}
                  onChange={(e) => setOnlyUserUploads(e.target.checked)}
                  disabled={discoveryMutation.isPending || compoundItems.length === 0}
                />
              }
              label="Only use my uploads"
              title="Skip the shared database entirely and align only against your own uploaded compounds and gene clusters."
            />

            <Button
              variant="contained"
              disabled={!canQuery}
              onClick={() => {
                setQueryOriginSmiles(selectedItem?.kind === "compound" ? selectedItem.smiles : null);
                discoveryMutation.mutate();
              }}
              startIcon={discoveryMutation.isPending ? <CircularProgress size={16} color="inherit" /> : undefined}
            >
              {discoveryMutation.isPending ? "Querying..." : "Query"}
            </Button>
          </Stack>
        </CardContent>
      </Card>

      {discoveryMutation.data && (
        <Card variant="outlined">
          <CardContent>
            <Typography component="h1" variant="subtitle1">
              Results
            </Typography>
            <Typography variant="body2" color="text.secondary" sx={{ mb: 1.5 }}>
              {discoveryMutation.data.candidatesConsidered} candidate(s) considered
              {discoveryMutation.data.candidatesSkipped > 0
                ? `, ${discoveryMutation.data.candidatesSkipped} skipped (unrecognized tokens)`
                : ""}
              . Showing top {discoveryMutation.data.results.length}
              {discoveryMutation.data.scoreMode
                ? ` (${DISCOVERY_SCORE_MODE_OPTIONS.find((o) => o.value === discoveryMutation.data?.scoreMode)?.label.toLowerCase()})`
                : ""}
              .
            </Typography>

            {discoveryMutation.data.results.length === 0 ? (
              <Typography variant="body2" color="text.secondary">
                No results.
              </Typography>
            ) : (
              <>
                <Stack direction="row" spacing={1.5} alignItems="center" flexWrap="wrap" sx={{ mb: 1 }}>
                  <ToggleButtonGroup
                    size="small"
                    exclusive
                    value={resultsView}
                    onChange={(_, value) => value && setResultsView(value)}
                  >
                    <ToggleButton value="pairwise">Pairwise</ToggleButton>
                    <ToggleButton value="msa">Multiple sequence alignment</ToggleButton>
                    <ToggleButton value="compare">Compare compounds</ToggleButton>
                    <ToggleButton value="shape">Shape (PMI)</ToggleButton>
                  </ToggleButtonGroup>

                  <Typography variant="caption" color="text.secondary">
                    {selectedForMsa.size} of {discoveryMutation.data.results.length} selected
                  </Typography>
                  <Button size="small" onClick={() => setSelectedForMsa(new Set(discoveryMutation.data!.results.map((r) => r.entryId)))}>
                    Select all
                  </Button>
                  <Button size="small" onClick={() => setSelectedForMsa(new Set())}>
                    Clear
                  </Button>
                </Stack>

                {resultsView !== "compare" && resultsView !== "shape" && (
                  <Typography variant="caption" color="text.secondary" sx={{ display: "block", mb: 1.5 }}>
                    Cell shading reflects each unit's structural similarity to the query's aligned unit in that
                    column. Darker means a closer match, so a sequence's weak spots stand out at a glance.
                  </Typography>
                )}

                <Divider sx={{ mb: 1 }} />

                {resultsView === "pairwise" ? (
                  discoveryMutation.data.results.map((result, idx) => (
                    <ResultRow
                      key={result.entryId}
                      result={result}
                      rank={idx + 1}
                      selectedForMsa={selectedForMsa.has(result.entryId)}
                      onToggleSelectedForMsa={() => toggleSelectedForMsa(result.entryId)}
                    />
                  ))
                ) : resultsView === "msa" ? (
                  selectedResults.length === 0 ? (
                    <Typography variant="body2" color="text.secondary">
                      Select at least one result above to compute a multiple sequence alignment.
                    </Typography>
                  ) : msaQuery.isLoading ? (
                    <CircularProgress size={20} />
                  ) : msaQuery.error ? (
                    <Alert severity="error">
                      {(msaQuery.error as Error).message || "Failed to compute the multiple sequence alignment."}
                    </Alert>
                  ) : (
                    <Box sx={{ py: 1 }}>
                      <AlignmentGrid rows={msaRows} />
                      <Box sx={{ mt: 1 }}>
                        <DownloadSvgButton rows={msaRows} filename="msa_alignment.svg" />
                      </Box>
                    </Box>
                  )
                ) : resultsView === "compare" ? (
                  selectedResults.length === 0 ? (
                    <Typography variant="body2" color="text.secondary">
                      Select at least one result above to compare compounds.
                    </Typography>
                  ) : (
                    <Box sx={{ py: 1 }}>
                      <WorkspaceCompare results={selectedResults} queryOriginSmiles={queryOriginSmiles} />
                    </Box>
                  )
                ) : selectedResults.length === 0 ? (
                  <Typography variant="body2" color="text.secondary">
                    Select at least one result above to analyze molecular shape.
                  </Typography>
                ) : (
                  <Box sx={{ py: 1 }}>
                    <WorkspaceShape results={selectedResults} queryOriginSmiles={queryOriginSmiles} />
                  </Box>
                )}
              </>
            )}
          </CardContent>
        </Card>
      )}
    </Box>
  );
};
