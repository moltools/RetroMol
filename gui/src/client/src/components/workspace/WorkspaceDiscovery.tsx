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
import FormControlLabel from "@mui/material/FormControlLabel";
import IconButton from "@mui/material/IconButton";
import ListSubheader from "@mui/material/ListSubheader";
import MenuItem from "@mui/material/MenuItem";
import TextField from "@mui/material/TextField";
import ToggleButton from "@mui/material/ToggleButton";
import ToggleButtonGroup from "@mui/material/ToggleButtonGroup";
import Tooltip from "@mui/material/Tooltip";
import DeleteIcon from "@mui/icons-material/Delete";
import UploadFileIcon from "@mui/icons-material/UploadFile";
import ViewIcon from "@mui/icons-material/Visibility";
import WarningAmberIcon from "@mui/icons-material/WarningAmber";
import { useQuery } from "@tanstack/react-query";
import { Session, SessionItem, DiscoveryQueryItem, DiscoveryQueryItemSchema } from "../../features/session/types";
import { deleteSessionItem } from "../../features/session/api";
import { reconstructCompound } from "../../features/reconstruction/api";
import type { PrimarySequenceItem } from "../../features/reconstruction/types";
import { reconstructGeneCluster } from "../../features/clusters/api";
import type { ClusterPrimarySequence } from "../../features/clusters/types";
import { submitDiscoveryQuery } from "../../features/discovery/api";
import {
  DISCOVERY_SCORE_MODE_OPTIONS,
  type DiscoveryEntryType,
  type DiscoveryScoreMode,
} from "../../features/discovery/types";
import { useNotifications } from "../NotificationProvider";
import { MotifName } from "../MotifName";
import { horizontalScrollSx } from "../../theme/scrollbarSx";
import { SequenceEditor, type SequenceBlock } from "./SequenceEditor";
import { DialogViewDiscoveryQuery } from "./DialogViewDiscoveryQuery";

type WorkspaceDiscoveryProps = {
  session: Session;
  setSession: React.Dispatch<React.SetStateAction<Session | null>>;
};

const MAX_TOP_X = 200;
export const MAX_DISCOVERY_QUERY_ITEMS = 5;

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

// One row in the "Saved queries" list -- mirrors WorkspaceItemCard's status-chip
// language (queued/processing/done/error) so a query reads the same way an uploaded
// compound/cluster does, just without the score gauge or expand/collapse panel.
function DiscoveryQueryListItem({
  item,
  deleting,
  onView,
  onDelete,
}: {
  item: DiscoveryQueryItem;
  deleting: boolean;
  onView: () => void;
  onDelete: () => void;
}) {
  const isQueued = item.status === "queued";
  const isProcessing = item.status === "processing";
  const isDone = item.status === "done";
  const isError = item.status === "error";

  const flagChips = [
    item.flags.computeMsa ? "MSA" : null,
    item.flags.computeCompare ? "Compare" : null,
    item.flags.computePmi ? "PMI" : null,
  ].filter((label): label is string => label !== null);

  return (
    <Box
      onClick={onView}
      sx={(theme) => {
        const t = theme.vars || theme;
        return {
          borderRadius: 1,
          border: `1px solid ${t.palette.divider}`,
          p: 1.5,
          display: "flex",
          alignItems: "center",
          gap: 1.5,
          cursor: "pointer",
          "&:hover": { boxShadow: 4 },
        };
      }}
    >
      <Box sx={{ flex: 1, minWidth: 0 }}>
        <Typography variant="body2" fontWeight={500} noWrap>
          {item.name}
        </Typography>
        <Stack direction="row" spacing={0.5} flexWrap="wrap" useFlexGap sx={{ mt: 0.5 }}>
          <Chip label={item.settings.entryType} size="small" sx={{ fontSize: "0.7rem", height: 18 }} />
          {flagChips.map((label) => (
            <Chip key={label} label={label} size="small" variant="outlined" sx={{ fontSize: "0.7rem", height: 18 }} />
          ))}
        </Stack>
      </Box>

      {isQueued && <Chip label="Queued" color="warning" size="small" sx={{ fontSize: "0.7rem", height: 20 }} />}
      {isProcessing && <CircularProgress size={16} thickness={4} />}
      {isDone && <Chip label="Ready" color="success" size="small" sx={{ fontSize: "0.7rem", height: 20 }} />}
      {isError && (
        <Tooltip title={item.errorMessage || "An unknown error occurred."} placement="left" arrow>
          <Chip label="Error" color="error" size="small" sx={{ fontSize: "0.7rem", height: 20 }} />
        </Tooltip>
      )}

      <IconButton
        size="small"
        disabled={deleting}
        onClick={(e) => {
          e.stopPropagation();
          onView();
        }}
      >
        <ViewIcon fontSize="small" />
      </IconButton>
      <IconButton
        size="small"
        disabled={deleting}
        onClick={(e) => {
          e.stopPropagation();
          onDelete();
        }}
      >
        {deleting ? <CircularProgress size={16} thickness={4} /> : <DeleteIcon fontSize="small" />}
      </IconButton>
    </Box>
  );
}

export const WorkspaceDiscovery: React.FC<WorkspaceDiscoveryProps> = ({ session, setSession }) => {
  const { pushNotification } = useNotifications();

  const compoundItems = session.items.filter((item): item is SessionItem & { kind: "compound" } => item.kind === "compound");
  const clusterItems = session.items.filter((item): item is SessionItem & { kind: "cluster" } => item.kind === "cluster");
  const queryItems = session.items.filter((item): item is DiscoveryQueryItem => item.kind === "discoveryQuery");

  const [selectedItemId, setSelectedItemId] = React.useState<string>("");
  const selectedItem = session.items.find((item) => item.id === selectedItemId);
  const [blocks, setBlocks] = React.useState<SequenceBlock[]>([]);

  const [entryType, setEntryType] = React.useState<DiscoveryEntryType>("compound");
  const [scoreMode, setScoreMode] = React.useState<DiscoveryScoreMode>("subsequence");
  const [n, setN] = React.useState<number>(100);
  const [topX, setTopX] = React.useState<number>(20);
  const [includeUserUploads, setIncludeUserUploads] = React.useState<boolean>(false);
  const [onlyUserUploads, setOnlyUserUploads] = React.useState<boolean>(false);
  const [computeMsa, setComputeMsa] = React.useState<boolean>(false);
  const [computeCompare, setComputeCompare] = React.useState<boolean>(false);
  const [computePmi, setComputePmi] = React.useState<boolean>(false);
  const [submitting, setSubmitting] = React.useState<boolean>(false);

  const [deletingIds, setDeletingIds] = React.useState<Set<string>>(new Set());
  const [viewingItemId, setViewingItemId] = React.useState<string | null>(null);
  const [uploadedViewItem, setUploadedViewItem] = React.useState<DiscoveryQueryItem | null>(null);

  // Clean up deletingIds when session items change (mirrors WorkspaceUpload).
  React.useEffect(() => {
    setDeletingIds((prev) => {
      const liveIds = new Set(session.items.map((it) => it.id));
      const next = new Set<string>();
      prev.forEach((id) => {
        if (liveIds.has(id)) next.add(id);
      });
      return next;
    });
  }, [session.items]);

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

  const handlePickNames = (names: string[]) => {
    setBlocks(blocksFromNames(names));
  };

  const maxTopX = Math.max(1, Math.min(n, MAX_TOP_X));
  const atQueryCap = queryItems.length >= MAX_DISCOVERY_QUERY_ITEMS;
  const canQuery = blocks.length > 0 && !submitting && !atQueryCap;

  const queryOriginSmiles = selectedItem?.kind === "compound" ? selectedItem.smiles : null;

  const handleSubmitQuery = async () => {
    if (!canQuery) return;
    setSubmitting(true);

    const preview = blocks.slice(0, 4).map((b) => b.name).join(" ");
    const name = blocks.length > 4 ? `${preview} …` : preview || "Discovery query";

    try {
      const item = await submitDiscoveryQuery({
        sessionId: session.sessionId,
        name,
        primarySequence: blocks.map((b) => b.name),
        entryType,
        scoreMode,
        n,
        topX,
        includeUserUploads: includeUserUploads || onlyUserUploads,
        onlyUserUploads,
        queryOriginSmiles,
        flags: { computeMsa, computeCompare, computePmi },
      });
      setSession((prev) => (prev ? { ...prev, items: [...prev.items, item] } : prev));
      setViewingItemId(item.id);
    } catch (err) {
      const msg = err instanceof Error ? err.message : String(err);
      pushNotification(`Failed to submit discovery query: ${msg}`, "error");
    } finally {
      setSubmitting(false);
    }
  };

  const handleDeleteQueryItem = async (id: string) => {
    setDeletingIds((prev) => new Set(prev).add(id));
    try {
      await deleteSessionItem(session.sessionId, id);
      // SSE refresh will remove the item from session; nothing else to do here.
    } catch (err) {
      const msg = err instanceof Error ? err.message : String(err);
      pushNotification(`Failed to delete query: ${msg}`, "error");
      setDeletingIds((prev) => {
        const n = new Set(prev);
        n.delete(id);
        return n;
      });
    }
  };

  const handleUploadResultFile = async (event: React.ChangeEvent<HTMLInputElement>) => {
    const file = event.target.files?.[0];
    event.target.value = ""; // allow re-selecting the same file later
    if (!file) return;

    try {
      const text = await file.text();
      const parsed = DiscoveryQueryItemSchema.parse(JSON.parse(text));
      setUploadedViewItem(parsed);
      setViewingItemId(null);
    } catch (err) {
      const msg = err instanceof Error ? err.message : String(err);
      pushNotification(`Failed to load result file: ${msg}`, "error");
    }
  };

  const viewingSessionItem = viewingItemId ? queryItems.find((it) => it.id === viewingItemId) ?? null : null;
  const activeViewItem = uploadedViewItem ?? viewingSessionItem;

  const handleCloseViewDialog = () => {
    setViewingItemId(null);
    setUploadedViewItem(null);
  };

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

          <SequenceEditor blocks={blocks} onChange={setBlocks} disabled={submitting} />
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
              disabled={submitting}
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
              disabled={submitting}
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
              disabled={submitting}
            />

            <TextField
              size="small"
              type="number"
              label="Show top X"
              value={topX}
              onChange={(e) => setTopX(Math.max(1, Math.min(maxTopX, Number(e.target.value) || 1)))}
              slotProps={{ htmlInput: { min: 1, max: maxTopX } }}
              sx={{ width: 140 }}
              disabled={submitting}
            />

            <FormControlLabel
              control={
                <Checkbox
                  size="small"
                  checked={includeUserUploads || onlyUserUploads}
                  onChange={(e) => setIncludeUserUploads(e.target.checked)}
                  disabled={submitting || compoundItems.length === 0 || onlyUserUploads}
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
                  disabled={submitting || compoundItems.length === 0}
                />
              }
              label="Only use my uploads"
              title="Skip the shared database entirely and align only against your own uploaded compounds and gene clusters."
            />
          </Stack>

          <Typography variant="subtitle2" sx={{ mt: 2.5, mb: 0.5 }}>
            Precompute for this query
          </Typography>
          <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
            Pick which extra views to compute up front so they open instantly once the query is done. A view you
            didn't flag here is disabled in the results, not silently missing -- run another query with it enabled
            if you need it later.
          </Typography>

          <Stack direction="row" spacing={2} alignItems="center" flexWrap="wrap">
            <FormControlLabel
              control={
                <Checkbox size="small" checked={computeMsa} onChange={(e) => setComputeMsa(e.target.checked)} disabled={submitting} />
              }
              label="Compute MSA"
            />
            <FormControlLabel
              control={
                <Checkbox
                  size="small"
                  checked={computeCompare}
                  onChange={(e) => setComputeCompare(e.target.checked)}
                  disabled={submitting}
                />
              }
              label="Compute compound comparison"
            />
            <FormControlLabel
              control={
                <Checkbox size="small" checked={computePmi} onChange={(e) => setComputePmi(e.target.checked)} disabled={submitting} />
              }
              label={
                <Stack direction="row" spacing={0.5} alignItems="center">
                  <span>Compute shape (PMI)</span>
                  <Tooltip title="Conformer search is slow -- this can take noticeably longer than the rest of the query." arrow>
                    <WarningAmberIcon fontSize="inherit" color="warning" />
                  </Tooltip>
                </Stack>
              }
            />
          </Stack>
          {computePmi && (
            <Alert severity="warning" sx={{ mt: 1 }}>
              Shape (PMI) runs a conformer search per compound and can take noticeably longer than the rest of the
              query. It's routed to a separate queue so it won't block other users' faster queries, but this query
              itself will stay "processing" until it's done.
            </Alert>
          )}

          <Stack direction="row" spacing={2} alignItems="center" sx={{ mt: 2.5 }}>
            <Button
              variant="contained"
              disabled={!canQuery}
              onClick={handleSubmitQuery}
              startIcon={submitting ? <CircularProgress size={16} color="inherit" /> : undefined}
            >
              {submitting ? "Submitting…" : "Query"}
            </Button>
            {atQueryCap && (
              <Typography variant="caption" color="text.secondary">
                You've reached the {MAX_DISCOVERY_QUERY_ITEMS} saved query limit -- delete one below to run another.
              </Typography>
            )}
          </Stack>
        </CardContent>
      </Card>

      <Card variant="outlined">
        <CardContent>
          <Stack direction="row" alignItems="center" justifyContent="space-between" flexWrap="wrap" gap={1} sx={{ mb: 1.5 }}>
            <Typography component="h1" variant="subtitle1">
              Saved queries ({queryItems.length}/{MAX_DISCOVERY_QUERY_ITEMS})
            </Typography>
            <Button size="small" variant="text" component="label" startIcon={<UploadFileIcon fontSize="small" />}>
              Load result file
              <input type="file" accept="application/json" hidden onChange={handleUploadResultFile} />
            </Button>
          </Stack>
          <Typography variant="caption" color="text.secondary" sx={{ display: "block", mb: 1.5 }}>
            Queries stay here (and survive switching tabs or reloading) until you delete them, up to{" "}
            {MAX_DISCOVERY_QUERY_ITEMS} at a time. "Load result file" opens a previously downloaded result for
            viewing only -- it doesn't count against this limit.
          </Typography>

          {queryItems.length === 0 ? (
            <Typography variant="body2" color="text.secondary">
              No saved queries yet. Configure and run a query above.
            </Typography>
          ) : (
            <Stack spacing={1}>
              {queryItems.map((item) => (
                <DiscoveryQueryListItem
                  key={item.id}
                  item={item}
                  deleting={deletingIds.has(item.id)}
                  onView={() => {
                    setUploadedViewItem(null);
                    setViewingItemId(item.id);
                  }}
                  onDelete={() => handleDeleteQueryItem(item.id)}
                />
              ))}
            </Stack>
          )}
        </CardContent>
      </Card>

      {activeViewItem && (
        <DialogViewDiscoveryQuery
          item={activeViewItem}
          sessionId={uploadedViewItem ? null : session.sessionId}
          open={!!activeViewItem}
          onClose={handleCloseViewDialog}
        />
      )}
    </Box>
  );
};
