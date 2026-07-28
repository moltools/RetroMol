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
import Tooltip from "@mui/material/Tooltip";
import DeleteIcon from "@mui/icons-material/Delete";
import ViewIcon from "@mui/icons-material/Visibility";
import ExpandMoreIcon from "@mui/icons-material/ExpandMore";
import ExpandLessIcon from "@mui/icons-material/ExpandLess";
import CircularProgress from "@mui/material/CircularProgress";
import { useQuery } from "@tanstack/react-query";
import { Gauge } from "@mui/x-charts/Gauge";
import { Session, SessionItem } from "../../features/session/types";
import { alpha } from "@mui/material/styles";
import type { Theme } from "@mui/material/styles";
import { DialogViewItem } from "./DialogViewItem";
import { PrimarySequenceRows, usePrimarySequenceEditor } from "./PrimarySequenceEditor";
import { ClusterReadoutRows } from "./ClusterReadoutRows";
import { reconstructCompound } from "../../features/reconstruction/api";
import { getClusterReadout } from "../../features/clusters/api";
import { useTick } from "../../hooks/useTick";

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

  const [openViewItem, setOpenViewItem] = React.useState(false);
  const [expanded, setExpanded] = React.useState(false);
  const [selectedTags, setSelectedTags] = React.useState<number[]>([]);

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

  const handleToggleMotif = (tags: number[]) => {
    setSelectedTags((prev) => {
      const allSelected = tags.every((tag) => prev.includes(tag));
      if (allSelected) return prev.filter((tag) => !tags.includes(tag));
      return Array.from(new Set([...prev, ...tags]));
    });
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

            <Stack direction="column" spacing={0.5} sx={{ minWidth: 0, flex: 1 }}>
              <Stack direction="row" spacing={0.5} alignItems="center" sx={{ minWidth: 0 }}>
                <Stack direction="row" alignItems="center" spacing={1} sx={{ minWidth: 0, flex: 1 }}>
                    <Typography
                      variant="body2"
                      fontWeight={500}
                      noWrap
                      sx={{
                        flex: 1,
                        minWidth: 0,
                        overflow: "hidden",
                        textOverflow: "ellipsis",
                        whiteSpace: "nowrap",
                      }}
                    >
                      {item.name}
                    </Typography>
                </Stack>
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
                      selectedTags={selectedTags}
                      onToggleMotif={handleToggleMotif}
                      labelWidth={110}
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
                      <Button size="small" variant="outlined" onClick={() => editor.setEditing(true)}>
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