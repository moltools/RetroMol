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
import IconButton from "@mui/material/IconButton";
import ListSubheader from "@mui/material/ListSubheader";
import MenuItem from "@mui/material/MenuItem";
import TextField from "@mui/material/TextField";
import ToggleButton from "@mui/material/ToggleButton";
import ToggleButtonGroup from "@mui/material/ToggleButtonGroup";
import Tooltip from "@mui/material/Tooltip";
import CheckIcon from "@mui/icons-material/Check";
import CloseIcon from "@mui/icons-material/Close";
import DeleteIcon from "@mui/icons-material/Delete";
import EditIcon from "@mui/icons-material/Edit";
import SettingsIcon from "@mui/icons-material/Settings";
import UploadFileIcon from "@mui/icons-material/UploadFile";
import ViewIcon from "@mui/icons-material/Visibility";
import { alpha } from "@mui/material/styles";
import { useQuery } from "@tanstack/react-query";
import { Session, SessionItem, DiscoveryQueryItem, DiscoveryQueryItemSchema } from "../../features/session/types";
import { deleteSessionItem, renameSessionItem } from "../../features/session/api";
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
import { PrimarySequenceOverview } from "./PrimarySequenceEditor";
import { DialogViewDiscoveryQuery } from "./DialogViewDiscoveryQuery";
import {MinimalIconButton} from "../MinimalIconButton";

type WorkspaceDiscoveryProps = {
  session: Session;
  setSession: React.Dispatch<React.SetStateAction<Session | null>>;
};

const MAX_TOP_X = 200;
export const MAX_DISCOVERY_QUERY_ITEMS = 5;

function blocksFromNames(names: string[]): SequenceBlock[] {
  return names.map((name) => ({ id: crypto.randomUUID(), name }));
}

// Default query name -- a block-name preview reads as noise once there are a few
// queries in the list, so name by submission time instead ("Query 2026-08-09 14:32:05").
function formatQueryTimestamp(date: Date): string {
  const pad = (n: number) => String(n).padStart(2, "0");
  return `${date.getFullYear()}-${pad(date.getMonth() + 1)}-${pad(date.getDate())} ${pad(date.getHours())}:${pad(date.getMinutes())}:${pad(date.getSeconds())}`;
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
    <Box sx={{ ...horizontalScrollSx }}>
      <Box
        sx={{
          display: "flex",
          flexWrap: "nowrap",
          gap: 0.5,
          alignItems: "center",
          position: "relative",
          width: "max-content",

          "&::before": {
            content: '""',
            position: "absolute",
            left: 0,
            right: 0,
            top: "50%",
            height: "2px",
            backgroundColor: "divider",
            zIndex: 0,
          },
        }}
      >
        {names.map((name, idx) => (
          <Box
            key={`${name}-${idx}`}
            sx={{
              px: 1,
              py: 0.5,
              borderRadius: 1,
              border: "1px solid",
              borderColor: "divider",
              bgcolor: "background.paper",
              color: "text.primary",
              fontSize: "0.8rem",
              whiteSpace: "nowrap",
              flexShrink: 0,
              zIndex: 1,
            }}
          >
            <MotifName name={name} />
          </Box>
        ))}
      </Box>
    </Box>
  );
}

// One row in the "Query results" list -- mirrors WorkspaceItemCard's look (selection
// checkbox, tinted background, inline rename) and status-chip language
// (queued/processing/done/error), just without the score gauge or expand/collapse panel.
function DiscoveryQueryListItem({
  session,
  setSession,
  item,
  selected,
  deleting,
  onToggleSelect,
  onView,
  onDelete,
}: {
  session: Session;
  setSession: React.Dispatch<React.SetStateAction<Session | null>>;
  item: DiscoveryQueryItem;
  selected: boolean;
  deleting: boolean;
  onToggleSelect: (id: string) => void;
  onView: () => void;
  onDelete: () => void;
}) {
  const { pushNotification } = useNotifications();

  const [editingName, setEditingName] = React.useState(false);
  const [nameDraft, setNameDraft] = React.useState(item.name);
  const [savingName, setSavingName] = React.useState(false);

  // Keep the draft in sync with the persisted name as long as we're not mid-edit.
  React.useEffect(() => {
    if (!editingName) setNameDraft(item.name);
  }, [item.name, editingName]);

  const isQueued = item.status === "queued";
  const isProcessing = item.status === "processing";
  const isDone = item.status === "done";
  const isError = item.status === "error";

  const flagChips = [
    item.flags.computeMsa ? "MSA" : null,
    item.flags.computeCompare ? "Compare" : null,
    item.flags.computeEnrichment ? "Enrichment" : null,
  ].filter((label): label is string => label !== null);

  const handleToggle = (e?: React.SyntheticEvent) => {
    if (e) e.stopPropagation();
    if (deleting) return;
    onToggleSelect(item.id);
  };

  const handleStartEditName = (e: React.SyntheticEvent) => {
    e.stopPropagation();
    if (deleting) return;
    setNameDraft(item.name);
    setEditingName(true);
  };

  const handleCancelEditName = () => {
    setNameDraft(item.name);
    setEditingName(false);
  };

  const handleSaveName = async () => {
    const trimmed = nameDraft.trim();
    if (!trimmed) {
      pushNotification("Name cannot be empty.", "error");
      return;
    }
    if (trimmed === item.name) {
      setEditingName(false);
      return;
    }

    setSavingName(true);
    try {
      const nextSession = await renameSessionItem(session, item.id, trimmed);
      setSession(() => nextSession);
      setEditingName(false);
    } catch (err) {
      const msg = err instanceof Error ? err.message : String(err);
      pushNotification(`Failed to rename query: ${msg}`, "error");
    } finally {
      setSavingName(false);
    }
  };

  return (
    <Stack
      onClick={handleToggle}
      direction="row"
      alignItems="center"
      sx={(theme) => {
        const t = theme.vars || theme;
        return {
          borderRadius: 1,
          border: `1px solid ${selected ? t.palette.primary.main : "transparent"}`,
          p: 1.5,
          gap: 1.5,
          // Matches WorkspaceItemCard's rendered height (driven there by its 70px score
          // gauge), so the two lists read as one consistent row height across tabs.
          minHeight: 96,
          cursor: "pointer",
          "&:hover": { boxShadow: 10 },
          backgroundColor: selected ? alpha("#000000", 0.04) : alpha("#000000", 0.02),
          ...theme.applyStyles("dark", { backgroundColor: selected ? alpha("#ffffff", 0.06) : alpha("#ffffff", 0.03) }),
        };
      }}
    >
      <Checkbox
        size="small"
        checked={selected}
        disabled={deleting}
        onClick={(e) => {
          e.stopPropagation();
          onToggleSelect(item.id);
        }}
      />

      <Box sx={{ flex: 1, minWidth: 0 }}>
        {editingName ? (
          <Stack
            direction="row"
            alignItems="center"
            spacing={0.5}
            sx={{ minWidth: 0, width: "fit-content", maxWidth: "100%" }}
            onClick={(e) => e.stopPropagation()}
          >
            <TextField
              value={nameDraft}
              onChange={(e) => setNameDraft(e.target.value)}
              onKeyDown={(e) => {
                if (e.key === "Enter") {
                  e.preventDefault();
                  handleSaveName();
                } else if (e.key === "Escape") {
                  e.preventDefault();
                  handleCancelEditName();
                }
              }}
              size="small"
              variant="standard"
              autoFocus
              disabled={savingName}
              sx={{ minWidth: 0, width: "auto" }}
            />
            <Tooltip title="Save name" arrow>
              <span>
                <IconButton
                  size="small"
                  onClick={handleSaveName}
                  disabled={savingName || !nameDraft.trim()}
                  sx={{ p: 0.25, flexShrink: 0 }}
                >
                  {savingName ? <CircularProgress size={14} /> : <CheckIcon fontSize="small" />}
                </IconButton>
              </span>
            </Tooltip>
            <Tooltip title="Cancel" arrow>
              <span>
                <IconButton size="small" onClick={handleCancelEditName} disabled={savingName} sx={{ p: 0.25, flexShrink: 0 }}>
                  <CloseIcon fontSize="small" />
                </IconButton>
              </span>
            </Tooltip>
          </Stack>
        ) : (
          <Stack direction="row" alignItems="center" spacing={0.5} sx={{ minWidth: 0 }}>
            <Typography
              variant="body2"
              fontWeight={500}
              noWrap
              sx={{ minWidth: 0, overflow: "hidden", textOverflow: "ellipsis", whiteSpace: "nowrap" }}
            >
              {item.name}
            </Typography>
            <Tooltip title="Rename" arrow>
              <span>
                <MinimalIconButton onClick={handleStartEditName} disabled={deleting}>
                  <EditIcon sx={{ fontSize: "0.9rem", transform: "translateY(-2px)" }} />
                </MinimalIconButton>
              </span>
            </Tooltip>
          </Stack>
        )}

        <Stack direction="row" spacing={0.5} flexWrap="wrap" useFlexGap sx={{ mt: 0.5 }}>
          <Chip label={item.settings.entryType} size="small" sx={{ fontSize: "0.7rem", height: 18 }} />
          {flagChips.map((label) => (
            <Chip key={label} label={label} size="small" variant="outlined" sx={{ fontSize: "0.7rem", height: 18 }} />
          ))}
        </Stack>
      </Box>

      {deleting && (
        <>
          <Chip label="Deleting..." size="small" sx={{ fontSize: "0.7rem", height: 20 }} />
          <CircularProgress size={16} thickness={4} />
        </>
      )}
      {isQueued && <Chip label="Queued" color="warning" size="small" sx={{ fontSize: "0.7rem", height: 20 }} />}
      {isProcessing && <CircularProgress size={16} thickness={4} />}
      {isDone && <Chip label="Ready" color="success" size="small" sx={{ fontSize: "0.7rem", height: 20 }} />}
      {isError && (
        <Tooltip title={item.errorMessage || "An unknown error occurred."} placement="left" arrow>
          <Chip label="Error" color="error" size="small" sx={{ fontSize: "0.7rem", height: 20 }} />
        </Tooltip>
      )}

      <Tooltip title="View details" arrow>
        <span>
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
        </span>
      </Tooltip>
      <Tooltip title="Delete" arrow>
        <span>
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
        </span>
      </Tooltip>
    </Stack>
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
  const [scoreMode, setScoreMode] = React.useState<DiscoveryScoreMode>("longest_sequence");
  const [n, setN] = React.useState<number>(1000);
  const [topX, setTopX] = React.useState<number>(20);
  const [includeUserUploads, setIncludeUserUploads] = React.useState<boolean>(false);
  const [onlyUserUploads, setOnlyUserUploads] = React.useState<boolean>(false);
  const [computeMsa, setComputeMsa] = React.useState<boolean>(true);
  const [computeCompare, setComputeCompare] = React.useState<boolean>(false);
  const [computeEnrichment, setComputeEnrichment] = React.useState<boolean>(true);
  const [submitting, setSubmitting] = React.useState<boolean>(false);
  const [queryOptionsOpen, setQueryOptionsOpen] = React.useState<boolean>(false);

  const [deletingIds, setDeletingIds] = React.useState<Set<string>>(new Set());
  const [selectedQueryIds, setSelectedQueryIds] = React.useState<Set<string>>(new Set());
  const [viewingItemId, setViewingItemId] = React.useState<string | null>(null);
  const [uploadedViewItem, setUploadedViewItem] = React.useState<DiscoveryQueryItem | null>(null);

  // Clean up deletingIds/selectedQueryIds when session items change (mirrors WorkspaceUpload).
  React.useEffect(() => {
    const liveIds = new Set(session.items.map((it) => it.id));
    const keepLive = (prev: Set<string>) => {
      const next = new Set<string>();
      prev.forEach((id) => {
        if (liveIds.has(id)) next.add(id);
      });
      return next;
    };
    setDeletingIds(keepLive);
    setSelectedQueryIds(keepLive);
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

  const anyQuerySelected = selectedQueryIds.size > 0;
  const allQueriesSelected = queryItems.length > 0 && selectedQueryIds.size === queryItems.length;

  const queryOriginSmiles = selectedItem?.kind === "compound" ? selectedItem.smiles : null;

  const handleSubmitQuery = async () => {
    if (!canQuery) return;
    setSubmitting(true);

    const name = `Query ${formatQueryTimestamp(new Date())}`;

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
        flags: { computeMsa, computeCompare, computeEnrichment },
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

  // Query-selection helpers -- mirrors WorkspaceUpload's select-all/clear/delete-selected trio.
  const toggleSelectQueryItem = (id: string) => {
    if (deletingIds.has(id)) return;
    setSelectedQueryIds((prev) => {
      const next = new Set(prev);
      if (next.has(id)) next.delete(id);
      else next.add(id);
      return next;
    });
  };

  const handleSelectAllQueries = () => {
    if (!queryItems.length) return;
    setSelectedQueryIds(new Set(queryItems.map((item) => item.id)));
  };

  const handleClearQuerySelection = () => {
    setSelectedQueryIds(new Set());
  };

  const handleDeleteSelectedQueries = async () => {
    if (selectedQueryIds.size === 0) return;

    const ids = Array.from(selectedQueryIds);
    setSelectedQueryIds(new Set());
    setDeletingIds((prev) => {
      const next = new Set(prev);
      ids.forEach((id) => next.add(id));
      return next;
    });

    try {
      for (const id of ids) {
        await deleteSessionItem(session.sessionId, id);
      }
      // SSE will update the session state
    } catch (err) {
      const msg = err instanceof Error ? err.message : String(err);
      pushNotification(`Failed to delete selected queries: ${msg}`, "error");
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
              {reconstructionQuery.data &&
                reconstructionQuery.data.dbMatchingSequences.length === 0 &&
                reconstructionQuery.data.reconstructions.length === 0 && (
                  <Typography variant="body2" color="text.secondary">
                    No reconstructed primary sequences found for this compound.
                  </Typography>
                )}

              {reconstructionQuery.data && reconstructionQuery.data.dbMatchingSequences.length > 0 && (
                <>
                  <Typography variant="caption" color="text.secondary" sx={{ display: "block", mb: 0.5 }}>
                    Database-matching sequence -- query with this for a fingerprint
                    guaranteed comparable to what's actually stored. Nothing is
                    filtered out of it (tailoring events like glycosylation are
                    included), so a compound with more than one biosynthetic chain
                    can also be queried with just one of them below.
                  </Typography>
                  <Stack spacing={1.5} sx={{ mb: 2 }}>
                    <PrimarySequenceOverview
                      data={reconstructionQuery.data.dbMatchingSequences}
                      renderAction={(sequence) => (
                        <Button
                          size="small"
                          variant="outlined"
                          onClick={() => handlePickNames(sequence.map(([name]) => name))}
                        >
                          Use this
                        </Button>
                      )}
                    />
                  </Stack>
                </>
              )}

              <Stack spacing={1}>
                {(reconstructionQuery.data?.reconstructions ?? []).map((reconstruction, idx) => {
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
                        minWidth: 0,
                      }}
                    >
                      <Typography variant="caption" color="text.secondary" sx={{ flexShrink: 0 }}>
                        {`reconstructed backbone ${idx+1}`}
                      </Typography>
                      {override && (
                        <Chip label="Edited" size="small" color="info" variant="outlined" sx={{ fontSize: "0.7rem", flexShrink: 0 }} />
                      )}
                      <Box sx={{ flex: 1, minWidth: 0, transform: "translateY(5px)" }}>
                        <ReconstructionPreview sequence={effectiveSequence} />
                      </Box>
                      <Button
                        size="small"
                        variant="outlined"
                        onClick={() => handlePickNames(effectiveSequence.map(([name]) => name))}
                        sx={{ flexShrink: 0 }}
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
                      minWidth: 0,
                    }}
                  >
                    <Typography variant="caption" color="text.secondary" sx={{ flexShrink: 0 }}>
                      {region.id}
                    </Typography>
                    <Box sx={{ flex: 1, minWidth: 0, transform: "translateY(5px)" }}>
                      <NamesPreview names={region.primary_sequence} />
                    </Box>
                    <Button
                      size="small"
                      variant="outlined"
                      onClick={() => handlePickNames(region.primary_sequence)}
                      sx={{ flexShrink: 0 }}
                    >
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
          <Typography variant="body2" color="text.secondary" sx={{ mb: 1.5 }}>
            Choose to query either compounds, BGCs, or both. Additional query settings can be revealed using the cog wheel.
          </Typography>

          <Stack direction="row" spacing={1} alignItems="center" sx={{ mt: 1.5 }}>
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

            <Tooltip title={queryOptionsOpen ? "Hide query options" : "Show query options"} arrow>
              <span>
                <MinimalIconButton
                  size="medium"
                  onClick={() => setQueryOptionsOpen((prev) => !prev)}
                  disabled={submitting}
                  aria-expanded={queryOptionsOpen}
                  aria-label="Toggle query options"
                >
                  <SettingsIcon
                    fontSize="medium"
                    sx={{
                      transition: "transform 300ms ease",
                      transform: queryOptionsOpen ? "rotate(180deg)" : "rotate(0deg)",
                    }}
                  />
                </MinimalIconButton>
              </span>
            </Tooltip>
          </Stack>

          <Collapse in={queryOptionsOpen} unmountOnExit>
            <Stack spacing={2.5} sx={{ mt: 4 }}>
              <Stack direction="row" spacing={2} alignItems="flex-start" flexWrap="wrap">
                <TextField
                  select
                  size="small"
                  label="Score mode"
                  value={scoreMode}
                  onChange={(e) => setScoreMode(e.target.value as DiscoveryScoreMode)}
                  disabled={submitting}
                  helperText={
                    scoreMode === "subsequence"
                      ? "Scores by the query's own length, surfacing subsequence matches."
                      : "Normalizes by the longer sequence, favoring similarly-sized matches."
                  }
                  sx={{ width: 240 }}
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
                  helperText="Nearest neighbors pulled from the database before ranking."
                  sx={{ width: 190 }}
                  disabled={submitting}
                />

                <TextField
                  size="small"
                  type="number"
                  label="Show top X"
                  value={topX}
                  onChange={(e) => setTopX(Math.max(1, Math.min(maxTopX, Number(e.target.value) || 1)))}
                  slotProps={{ htmlInput: { min: 1, max: maxTopX } }}
                  helperText="How many top-ranked results to display."
                  sx={{ width: 170 }}
                  disabled={submitting}
                />
              </Stack>

              <Stack direction="row" spacing={2} alignItems="center" flexWrap="wrap">
                <Tooltip
                  title="Uploaded compounds and gene clusters compete for a spot among the nearest N candidates, then follow the usual top-X ranking."
                  arrow
                >
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
                  />
                </Tooltip>

                <Tooltip
                  title="Skip the shared database entirely and align only against your own uploaded compounds and gene clusters."
                  arrow
                >
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
                  />
                </Tooltip>
              </Stack>

              <Divider />

              <Box>
                <Typography variant="subtitle2" sx={{ mb: 0.5 }}>
                  Compute for this query
                </Typography>
                <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
                  Pick which extra views to compute up front so they open instantly once the query is done. A view
                  you didn't flag here is disabled in the results, not silently missing. Run another query with it
                  enabled if you need it later.
                </Typography>

                <Stack direction="row" spacing={2} alignItems="center" flexWrap="wrap">
                  <FormControlLabel
                    control={
                      <Checkbox
                        size="small"
                        checked={computeMsa}
                        onChange={(e) => setComputeMsa(e.target.checked)}
                        disabled={submitting}
                      />
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
                  <Tooltip
                    title="Test whether the nearest neighbors are enriched for any annotation (phylogeny, chemical class, ...) vs. the rest of the database. Recalculates in the results view if you narrow down the selection."
                    arrow
                  >
                    <FormControlLabel
                      control={
                        <Checkbox
                          size="small"
                          checked={computeEnrichment}
                          onChange={(e) => setComputeEnrichment(e.target.checked)}
                          disabled={submitting}
                        />
                      }
                      label="Compute enrichment"
                    />
                  </Tooltip>
                </Stack>
              </Box>
            </Stack>
          </Collapse>

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
              Query results ({queryItems.length}/{MAX_DISCOVERY_QUERY_ITEMS})
            </Typography>
            <Stack direction="row" spacing={1} flexWrap="wrap">
              <Button
                size="small"
                variant="text"
                onClick={handleSelectAllQueries}
                disabled={queryItems.length === 0 || allQueriesSelected}
              >
                Select all
              </Button>
              <Button size="small" variant="text" onClick={handleClearQuerySelection} disabled={!anyQuerySelected}>
                Clear selection
              </Button>
              <Button
                size="small"
                variant="text"
                color="error"
                onClick={handleDeleteSelectedQueries}
                disabled={!anyQuerySelected}
              >
                Delete selected
              </Button>
              <Button size="small" variant="text" component="label" startIcon={<UploadFileIcon fontSize="small" />}>
                Load result file
                <input type="file" accept="application/json" hidden onChange={handleUploadResultFile} />
              </Button>
            </Stack>
          </Stack>
          <Typography variant="body2" color="text.secondary" sx={{ display: "block", mb: 1.5 }}>
            Queries stay here (and survive switching tabs or reloading) until you delete them, up to{" "}
            {MAX_DISCOVERY_QUERY_ITEMS} at a time. "Load result file" opens a previously downloaded result for
            viewing only. It doesn't count against this limit.
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
                  session={session}
                  setSession={setSession}
                  item={item}
                  selected={selectedQueryIds.has(item.id)}
                  deleting={deletingIds.has(item.id)}
                  onToggleSelect={toggleSelectQueryItem}
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
          session={session}
          setSession={setSession}
          open={!!activeViewItem}
          onClose={handleCloseViewDialog}
        />
      )}
    </Box>
  );
};
