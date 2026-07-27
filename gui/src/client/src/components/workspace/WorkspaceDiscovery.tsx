import React from "react";
import Box from "@mui/material/Box";
import Card from "@mui/material/Card";
import CardContent from "@mui/material/CardContent";
import Stack from "@mui/material/Stack";
import Typography from "@mui/material/Typography";
import Alert from "@mui/material/Alert";
import Button from "@mui/material/Button";
import Chip from "@mui/material/Chip";
import CircularProgress from "@mui/material/CircularProgress";
import Collapse from "@mui/material/Collapse";
import Divider from "@mui/material/Divider";
import MenuItem from "@mui/material/MenuItem";
import TextField from "@mui/material/TextField";
import ToggleButton from "@mui/material/ToggleButton";
import ToggleButtonGroup from "@mui/material/ToggleButtonGroup";
import MuiLink from "@mui/material/Link";
import IconButton from "@mui/material/IconButton";
import ExpandMoreIcon from "@mui/icons-material/ExpandMore";
import ExpandLessIcon from "@mui/icons-material/ExpandLess";
import { useQuery, useMutation } from "@tanstack/react-query";
import { useTheme } from "@mui/material/styles";
import { Session, SessionItem } from "../../features/session/types";
import { reconstructCompound } from "../../features/reconstruction/api";
import type { PrimarySequenceItem } from "../../features/reconstruction/types";
import { runDiscoveryQuery } from "../../features/discovery/api";
import {
  DISCOVERY_SCORE_MODE_OPTIONS,
  type DiscoveryEntryType,
  type DiscoveryResult,
  type DiscoveryScoreMode,
} from "../../features/discovery/types";
import { useNotifications } from "../NotificationProvider";
import { MotifName } from "../MotifName";
import { horizontalScrollSx } from "../../theme/scrollbarSx";
import { SequenceEditor, type SequenceBlock } from "./SequenceEditor";

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
  return (
    <Box sx={{ display: "flex", flexWrap: "nowrap", gap: 0.5, minWidth: 0, flex: "1 1 0%", ...horizontalScrollSx }}>
      {sequence.map(([name], idx) => (
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

function alignedCellSx(name: string | null, columnWidthCh: number) {
  return {
    px: 0.75,
    py: 0.4,
    borderRadius: 1,
    border: "1px dashed",
    borderColor: name === null ? "divider" : "primary.main",
    bgcolor: name === null ? "transparent" : "action.hover",
    fontSize: "0.75rem",
    width: `${columnWidthCh}ch`,
    boxSizing: "border-box" as const,
    textAlign: "center" as const,
    flexShrink: 0,
  };
}

// Renders the query/target rows together (rather than as two independent rows)
// so each column can be sized to fit its longer of the two names -- otherwise a
// short name like "X" lining up under a long one like "propionic acid" makes the
// alignment unreadable.
function AlignedPair({
  query,
  target,
  targetLabel,
}: {
  query: (string | null)[];
  target: (string | null)[];
  targetLabel: string;
}) {
  const length = Math.max(query.length, target.length);
  const columnWidths = Array.from({ length }, (_, idx) => {
    const q = query[idx]?.length ?? 1;
    const t = target[idx]?.length ?? 1;
    return Math.max(q, t, 1) + 1.5; // small buffer beyond the widest label
  });

  const renderPills = (sequence: (string | null)[]) => (
    <Stack direction="row" spacing={0.5} sx={{ flexWrap: "nowrap" }}>
      {columnWidths.map((width, idx) => {
        const name = sequence[idx] ?? null;
        return (
          <Box key={idx} sx={alignedCellSx(name, width)}>
            {name === null ? "–" : <MotifName name={name} />}
          </Box>
        );
      })}
    </Stack>
  );

  return (
    <Stack direction="row" spacing={0.5}>
      {/* Labels stay fixed outside the scroll area; both pill rows share one
          scroll container so they always scroll together, keeping columns aligned. */}
      <Stack spacing={0.5} sx={{ flexShrink: 0, justifyContent: "center" }}>
        <Typography variant="caption" sx={{ minWidth: 70, height: "1.5em", color: "text.secondary", fontWeight: 600 }}>
          Query
        </Typography>
        <Typography variant="caption" sx={{ minWidth: 70, height: "1.5em", color: "text.secondary", fontWeight: 600 }}>
          {targetLabel}
        </Typography>
      </Stack>
      <Box sx={{ minWidth: 0, flex: "1 1 0%", ...horizontalScrollSx }}>
        <Stack spacing={0.5}>
          {renderPills(query)}
          {renderPills(target)}
        </Stack>
      </Box>
    </Stack>
  );
}

function ResultRow({ result, rank }: { result: DiscoveryResult; rank: number }) {
  const [expanded, setExpanded] = React.useState(false);
  const theme = useTheme();

  return (
    <Box sx={{ borderBottom: "1px solid", borderColor: "divider", py: 1 }}>
      <Stack direction="row" alignItems="center" spacing={1.5}>
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

        <Chip label={result.type === "compound" ? "Compound" : "BGC"} size="small" sx={{ fontSize: "0.7rem" }} />

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
          <AlignedPair
            query={result.alignedQuery}
            target={result.alignedTarget}
            targetLabel={result.type === "compound" ? "Compound" : "BGC"}
          />
        </Box>
      </Collapse>
    </Box>
  );
}

export const WorkspaceDiscovery: React.FC<WorkspaceDiscoveryProps> = ({ session }) => {
  const { pushNotification } = useNotifications();

  const compoundItems = session.items.filter((item): item is SessionItem & { kind: "compound" } => item.kind === "compound");
  const clusterCount = session.items.length - compoundItems.length;

  const [selectedItemId, setSelectedItemId] = React.useState<string>("");
  const [blocks, setBlocks] = React.useState<SequenceBlock[]>([]);

  const [entryType, setEntryType] = React.useState<DiscoveryEntryType>("compound");
  const [scoreMode, setScoreMode] = React.useState<DiscoveryScoreMode>("subsequence");
  const [n, setN] = React.useState<number>(100);
  const [topX, setTopX] = React.useState<number>(20);

  const reconstructionQuery = useQuery({
    queryKey: ["reconstructCompound", session.sessionId, selectedItemId],
    queryFn: ({ signal }) => reconstructCompound(session.sessionId, selectedItemId, signal),
    enabled: selectedItemId.length > 0,
  });

  const discoveryMutation = useMutation({
    mutationFn: () =>
      runDiscoveryQuery({
        primarySequence: blocks.map((b) => b.name),
        entryType,
        scoreMode,
        n,
        topX,
      }),
    onError: (err) => {
      const msg = err instanceof Error ? err.message : String(err);
      pushNotification(`Discovery query failed: ${msg}`, "error");
    },
  });

  const handlePickReconstruction = (sequence: PrimarySequenceItem[]) => {
    setBlocks(blocksFromNames(sequence.map(([name]) => name)));
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
            Seed the editor below from one of your uploaded compounds, or build a sequence from scratch by adding blocks directly.
          </Typography>

          {clusterCount > 0 && (
            <Alert severity="info" sx={{ mb: 1.5 }}>
              {clusterCount} gene cluster upload(s) aren't shown here — primary-sequence picking isn't available for BGCs yet.
            </Alert>
          )}

          {compoundItems.length === 0 ? (
            <Typography variant="body2" color="text.secondary">
              No compounds uploaded yet. Import one from the Upload tab first.
            </Typography>
          ) : (
            <TextField
              select
              size="small"
              label="Compound"
              value={selectedItemId}
              onChange={(e) => setSelectedItemId(e.target.value)}
              sx={{ minWidth: 260 }}
            >
              {compoundItems.map((item) => (
                <MenuItem key={item.id} value={item.id} disabled={item.status !== "done"}>
                  {item.name} {item.status !== "done" ? `(${item.status})` : ""}
                </MenuItem>
              ))}
            </TextField>
          )}

          {selectedItemId && (
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
                {(reconstructionQuery.data ?? []).map((reconstruction, idx) => (
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
                    <ReconstructionPreview sequence={reconstruction.primary_sequence} />
                    <Button size="small" variant="outlined" onClick={() => handlePickReconstruction(reconstruction.primary_sequence)}>
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
            Drag blocks to reorder, click × to delete, or click + to add a new block at the end.
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

            <Button
              variant="contained"
              disabled={!canQuery}
              onClick={() => discoveryMutation.mutate()}
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
                <Divider sx={{ mb: 1 }} />
                {discoveryMutation.data.results.map((result, idx) => (
                  <ResultRow key={result.entryId} result={result} rank={idx + 1} />
                ))}
              </>
            )}
          </CardContent>
        </Card>
      )}
    </Box>
  );
};
