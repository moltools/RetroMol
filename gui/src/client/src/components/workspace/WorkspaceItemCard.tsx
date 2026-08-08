import React from "react";
import Alert from "@mui/material/Alert";
import Box from "@mui/material/Box";
import Button from "@mui/material/Button";
import Collapse from "@mui/material/Collapse";
import Divider from "@mui/material/Divider";
import Stack from "@mui/material/Stack";
import Typography from "@mui/material/Typography";
import Checkbox from "@mui/material/Checkbox";
import Chip from "@mui/material/Chip";
import IconButton from "@mui/material/IconButton";
import TextField from "@mui/material/TextField";
import Tooltip from "@mui/material/Tooltip";
import DeleteIcon from "@mui/icons-material/Delete";
import ViewIcon from "@mui/icons-material/Visibility";
import EditIcon from "@mui/icons-material/Edit";
import CheckIcon from "@mui/icons-material/Check";
import CloseIcon from "@mui/icons-material/Close";
import ExpandMoreIcon from "@mui/icons-material/ExpandMore";
import ExpandLessIcon from "@mui/icons-material/ExpandLess";
import CircularProgress from "@mui/material/CircularProgress";
import { useQuery } from "@tanstack/react-query";
import { Gauge } from "@mui/x-charts/Gauge";
import { Session, SessionItem } from "../../features/session/types";
import { renameSessionItem } from "../../features/session/api";
import { alpha } from "@mui/material/styles";
import type { Theme } from "@mui/material/styles";
import { DialogViewItem } from "./DialogViewItem";
import { PrimarySequenceRows, usePrimarySequenceEditor } from "./PrimarySequenceEditor";
import { ClusterReadoutRows } from "./ClusterReadoutRows";
import { reconstructCompound } from "../../features/reconstruction/api";
import { getClusterReadout } from "../../features/clusters/api";
import { useTick } from "../../hooks/useTick";
import { useNotifications } from "../NotificationProvider";
import { MinimalIconButton } from "../MinimalIconButton";

function getScoreColor(theme: Theme, value: number): string {
  const t = theme.vars || theme;
  if (value < 0.5) { return t.palette.error.main };
  if (value < 0.9) { return t.palette.warning.main };
  return t.palette.success.main;
};

type WorkspaceItemCardProps = {
  session: Session;
  setSession: React.Dispatch<React.SetStateAction<Session | null>>;
  item: SessionItem;
  selected: boolean;
  disabled?: boolean;
  onToggleSelect: (id: string) => void;
  onDelete: (id: string) => void;
};

// Helper to format "X ago"
function formatUpdatedAgo(updatedAt?: number): string {
  if (!updatedAt) return "Never updated";
  const now = Date.now();
  const diffMs = now - updatedAt;
  if (diffMs < 0) return "just now";

  const diffSec = Math.floor(diffMs / 1000);
  if (diffSec < 5) return "just now";
  if (diffSec < 60) return `${diffSec}s ago`;

  const diffMin = Math.floor(diffSec / 60);
  if (diffMin < 60) return `${diffMin}m ago`;

  const diffHours = Math.floor(diffMin / 60);
  if (diffHours < 24) return `${diffHours}h ago`;

  const diffDays = Math.floor(diffHours / 24);
  return `${diffDays}d ago`;
};

export const WorkspaceItemCard: React.FC<WorkspaceItemCardProps> = ({
  session,
  setSession,
  item,
  selected,
  disabled = false,
  onToggleSelect,
  onDelete,
}) => {
  const isCompound = item.kind === "compound"; // there are only two types: "compound" and "cluster"
  const itemScore = typeof item.score === "number" ? item.score : 0.0;

  const { pushNotification } = useNotifications();

  const [openViewItem, setOpenViewItem] = React.useState(false);
  const [expanded, setExpanded] = React.useState(false);

  const [editingName, setEditingName] = React.useState(false);
  const [nameDraft, setNameDraft] = React.useState(item.name);
  const [savingName, setSavingName] = React.useState(false);

  // Keep the draft in sync with the persisted name as long as we're not mid-edit
  // (e.g. another tab/session update coming through).
  React.useEffect(() => {
    if (!editingName) setNameDraft(item.name);
  }, [item.name, editingName]);

  // Re-render every 5s so "X ago" updates, via one shared timer for all cards
  useTick(5000);

  const isQueued = item.status === "queued";
  const showSpinner = item.status === "processing";
  const isError = item.status === "error";
  const isDone = item.status === "done";

  const handleToggle = (e?: React.SyntheticEvent) => {
    if (e) e.stopPropagation();
    if (disabled) return;
    onToggleSelect(item.id);
  };

  const handleOpenViewItem = (event: React.MouseEvent<HTMLButtonElement, MouseEvent>) => {
    event.currentTarget.blur(); // prevents 'Blocked aria-hidden on an element' warning
    setOpenViewItem(true);
  };

  const handleStartEditName = (e: React.SyntheticEvent) => {
    e.stopPropagation();
    if (disabled) return;
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
      pushNotification(`Failed to rename item: ${msg}`, "error");
    } finally {
      setSavingName(false);
    }
  };

  // Shares its query cache (and edit state machine) with DialogViewItem -- expanding
  // here and opening "View item" for the same compound don't refetch or diverge.
  const reconstructionQuery = useQuery({
    queryKey: ["reconstructCompound", session.sessionId, item.id],
    queryFn: ({ signal }) => reconstructCompound(session.sessionId, item.id, signal),
    enabled: expanded && isCompound,
  });
  const reconstructions = reconstructionQuery.data ?? null;
  const editor = usePrimarySequenceEditor(session, setSession, item, reconstructions);
  // See DialogViewItem's isUnordered -- an unordered "bag of motifs" result isn't a
  // sequence, so there's nothing meaningful to drag-and-drop reorder.
  const isUnordered = !!reconstructions?.length && reconstructions.every((r) => r.ordered === false);

  // Same idea as reconstructionQuery, but for gene clusters -- the parsed module
  // readout is computed server-side at submit time and only ever handed over via
  // this dedicated endpoint (item.payload is stripped from the session itself).
  const clusterReadoutQuery = useQuery({
    queryKey: ["getClusterReadout", session.sessionId, item.id],
    queryFn: ({ signal }) => getClusterReadout(session.sessionId, item.id, signal),
    enabled: expanded && !isCompound && isDone,
  });

  return (
    <>
      <Stack
        onClick={handleToggle}
        direction="column"
        sx={(theme) => {
          const t = theme.vars || theme;
          return {
            borderRadius: 1,
            border: `1px solid ${selected ? t.palette.primary.main : "transparent"}`,
            p: 1.5,
            display: "flex",
            gap: 1.5,
            cursor: "pointer",
            "&:hover": { boxShadow: 10 },
            backgroundColor: selected ? alpha("#000000", 0.04) : alpha("#000000", 0.02),
            ...theme.applyStyles("dark", { backgroundColor: selected ? alpha("#ffffff", 0.06) : alpha("#ffffff", 0.03) }),
          }
        }}
      >
        <Box
          sx={{
            display: "flex",
            alignItems: { xs: "flex-start", sm: "center" },
            justifyContent: { xs: "flex-start", sm: "space-between" },
            flexWrap: "wrap",
            gap: 1.5,
          }}
        >
          <Stack
            direction="row"
            spacing={1.5}
            alignItems="center"
            sx={{ flex: "1 1 260px", minWidth: 0 }}
          >
            <Checkbox
              size="small"
              checked={selected}
              disabled={disabled}
              onClick={(e) => {
                e.stopPropagation();
                onToggleSelect(item.id);
              }}
            />

            <Tooltip
              title={
                isCompound
                  ? "Coverage is the share of the molecule RetroMol could match to known building blocks. Low or 0% coverage means most (or all) of the structure falls outside its rule set, e.g. an unusual scaffold or a modification the current rules don't recognize. A low coverage leads to primary sequences that may be sparse, incomplete, or empty."
                  : "Coverage is the share of predicted modules RetroMol could confidently assign a substrate to. Low or 0% usually means PARAS couldn't confidently predict a substrate for the NRPS domains present, so those modules are marked unknown."
              }
              placement="right"
              arrow
            >
              <Gauge
                value={Math.round(itemScore * 100)}
                valueMin={0}
                valueMax={100}
                startAngle={-110}
                endAngle={110}
                width={70}
                height={70}
                innerRadius="70%"
                outerRadius="100%"
                sx={{
                  minWidth: 70,
                  cursor: "help",
                  "& text": {
                    fontSize: "0.65rem",
                    fontWeight: 600,
                  },
                  "& .MuiGauge-valueArc": {
                    fill: (theme) => getScoreColor(theme, itemScore),
                    transition: "stroke-dashoffset 0.3s ease",
                  },
                }}
                text={({ value }) => `${value}%`}
              />
            </Tooltip>

            <Stack direction="column" spacing={0.5} sx={{ minWidth: 0 }}>
              <Stack direction="row" spacing={0.5} alignItems="center" sx={{ minWidth: 0 }}>
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
                        <IconButton
                          size="small"
                          onClick={handleCancelEditName}
                          disabled={savingName}
                          sx={{ p: 0.25, flexShrink: 0 }}
                        >
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
                      sx={{
                        minWidth: 0,
                        overflow: "hidden",
                        textOverflow: "ellipsis",
                        whiteSpace: "nowrap",
                      }}
                    >
                      {item.name}
                    </Typography>
                    <Tooltip title="Rename" arrow>
                      <span>
                        <MinimalIconButton onClick={handleStartEditName} disabled={disabled}>
                          <EditIcon sx={{ fontSize: "0.9rem", transform: "translateY(-2px)" }} />
                        </MinimalIconButton>
                      </span>
                    </Tooltip>
                  </Stack>
                )}
              </Stack>

              <Typography variant="caption" color="text.secondary">
                Status updated {formatUpdatedAgo(item.updatedAt)}
              </Typography>
            </Stack>
          </Stack>

          <Stack
            direction="row"
            spacing={1}
            alignItems="center"
            useFlexGap
            sx={{
              flex: "0 1 auto",
              flexWrap: "wrap",
              justifyContent: { xs: "flex-start", sm: "flex-end" },
              maxWidth: "100%",
            }}
          >

            {disabled && (
              <>
                <Chip
                  label="Deleting..."
                  size="small"
                  sx={{ fontSize: "0.7rem", height: 20 }}
                />
                <CircularProgress size={16} thickness={4} />
              </>
            )}

            {isCompound ? (
              <Chip
                label="Compound"
                size="small"
                sx={{ fontSize: "0.7rem", height: 20 }}
              />
            ) : (
              <Chip
                label="BGC"
                size="small"
                sx={{ fontSize: "0.7rem", height: 20 }}
              />
            )}

            {isCompound && (
              <Chip
                label={item.matchStereochemistry ? "Stereo" : "Non-stereo"}
                size="small"
                sx={{ fontSize: "0.7rem", height: 20 }}
              />
            )}

            {item.kind === "cluster" && (
              <Chip
                label={`PARAS ≥ ${item.parasThreshold.toFixed(2)}`}
                size="small"
                title="Minimum PARAS prediction probability required to call an NRPS substrate"
                sx={{ fontSize: "0.7rem", height: 20 }}
              />
            )}

            {isQueued && (
              <Chip
                label="Queued"
                color="warning"
                size="small"
                sx={{ fontSize: "0.7rem", height: 20 }}
              />
            )}

            {showSpinner && (<CircularProgress size={16} thickness={4} />)}

            {isDone && (
              <Chip
                label="Ready"
                color="success"
                size="small"
                sx={{ fontSize: "0.7rem", height: 20 }}
              />
            )}

            {isError && (
              <Tooltip
                title={item.errorMessage || "An unknown error occurred."}
                placement="left"
                arrow
              >
                <Chip
                  label="Error"
                  color="error"
                  size="small"
                  sx={{ fontSize: "0.7rem", height: 20 }}
                />
              </Tooltip>
            )}
            <Tooltip title="View details" arrow>
              <IconButton
                size="small"
                disabled={disabled}
                onClick={(event) => {
                  event.stopPropagation();
                  if (disabled) return;
                  handleOpenViewItem(event);
                }}
              >
                <ViewIcon fontSize="small" />
              </IconButton>
            </Tooltip>
            <Tooltip title="Delete" arrow>
              <IconButton
                size="small"
                disabled={disabled}
                onClick={(e) => {
                  e.stopPropagation();
                  if (disabled) return;
                  onDelete(item.id);
                }}
              >
                <DeleteIcon fontSize="small" />
              </IconButton>
            </Tooltip>
            <Tooltip title={isCompound ? "Show primary sequences" : "Show parsed modules"} arrow>
              <span>
                <IconButton
                  size="small"
                  disabled={disabled}
                  onClick={(e) => {
                    e.stopPropagation();
                    if (disabled) return;
                    setExpanded((prev) => !prev);
                  }}
                >
                  {expanded ? <ExpandLessIcon fontSize="small" /> : <ExpandMoreIcon fontSize="small" />}
                </IconButton>
              </span>
            </Tooltip>
          </Stack>
        </Box>

        {isCompound && (
          <Collapse in={expanded} timeout="auto" unmountOnExit>
            <Box onClick={(e) => e.stopPropagation()} sx={{ pt: 0.5 }}>
              <Divider sx={{ mb: 1.5 }} />

              {reconstructionQuery.isLoading && <CircularProgress size={20} />}

              {reconstructionQuery.error && (
                <Alert severity="error">
                  {(reconstructionQuery.error as Error).message || "Failed to load reconstruction."}
                </Alert>
              )}

              {reconstructions && reconstructions.length === 0 && (
                <Typography variant="body2" color="text.secondary">
                  No reconstructed primary sequences found for this compound.
                </Typography>
              )}

              {reconstructions && reconstructions.length > 0 && (
                <Stack spacing={1.5}>
                  <Stack direction="column" spacing={1.5}>
                    <PrimarySequenceRows
                      item={item}
                      data={reconstructions}
                      state={editor}
                      selectedTags={[]}
                    />
                  </Stack>

                  <Stack direction="row" spacing={1}>
                    {editor.editing ? (
                      <>
                        <Button size="small" variant="text" color="inherit" onClick={editor.handleCancelEdit}>
                          Cancel
                        </Button>
                        <Button
                          size="small"
                          variant="contained"
                          onClick={editor.handleSaveAll}
                          disabled={!editor.anyDirty || editor.saving}
                          startIcon={editor.saving ? <CircularProgress size={14} color="inherit" /> : undefined}
                        >
                          {editor.saving ? "Saving..." : "Save changes"}
                        </Button>
                      </>
                    ) : (
                      <Button
                        size="small"
                        variant="contained"
                        onClick={() => editor.setEditing(true)}
                        disabled={isUnordered}
                      >
                        Edit sequences
                      </Button>
                    )}
                  </Stack>
                </Stack>
              )}
            </Box>
          </Collapse>
        )}

        {!isCompound && (
          <Collapse in={expanded} timeout="auto" unmountOnExit>
            <Box onClick={(e) => e.stopPropagation()} sx={{ pt: 0.5 }}>
              <Divider sx={{ mb: 1.5 }} />

              {isQueued && (
                <Typography variant="body2" color="text.secondary">
                  Waiting to be parsed...
                </Typography>
              )}

              {showSpinner && <CircularProgress size={20} />}

              {isError && (
                <Typography variant="body2" color="text.secondary">
                  Parsing failed -- see the error above for details.
                </Typography>
              )}

              {isDone && clusterReadoutQuery.isLoading && <CircularProgress size={20} />}

              {isDone && clusterReadoutQuery.error && (
                <Alert severity="error">
                  {(clusterReadoutQuery.error as Error).message || "Failed to load parsed gene cluster."}
                </Alert>
              )}

              {isDone && clusterReadoutQuery.data && (
                <ClusterReadoutRows readouts={clusterReadoutQuery.data.readouts} />
              )}
            </Box>
          </Collapse>
        )}
      </Stack>

      <DialogViewItem
        session={session}
        setSession={setSession}
        item={item}
        open={openViewItem}
        onClose={() => setOpenViewItem(false)}
      />
    </>
  );
};