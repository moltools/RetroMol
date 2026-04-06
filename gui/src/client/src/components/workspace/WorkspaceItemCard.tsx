import React from "react";
import Box from "@mui/material/Box";
import Stack from "@mui/material/Stack";
import Typography from "@mui/material/Typography";
import Checkbox from "@mui/material/Checkbox";
import Chip from "@mui/material/Chip";
import IconButton from "@mui/material/IconButton";
import TextField from "@mui/material/TextField";
import Tooltip from "@mui/material/Tooltip";
import ScienceIcon from "@mui/icons-material/Science";
import BiotechIcon from "@mui/icons-material/Biotech";
import DeleteIcon from "@mui/icons-material/Delete";
import CircularProgress from "@mui/material/CircularProgress";
import EditIcon from "@mui/icons-material/Edit";
import VisibilityIcon from "@mui/icons-material/Visibility";
import { SessionItem } from "../../features/session/types";
import { alpha } from "@mui/material/styles";
import type { Theme } from "@mui/material/styles";
import { ScoreBar } from "../ScoreBar";
import { ItemKindChip } from "../ItemKindChip";

type WorkspaceItemCardProps = {
  item: SessionItem;
  selected: boolean;
  onToggleSelect: (id: string) => void;
  onView: (id: string) => void;
  onDelete: (id: string) => void;
  onRename: (id: string, newName: string) => void;
}

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
}

function getScoreColor(theme: Theme, value: number): string {
  const t = theme.vars || theme;
  if (value < 0.5) { return t.palette.error.main };
  if (value < 0.9) { return t.palette.warning.main };
  return t.palette.success.main;
}

function getScoreTooltip(value: number): string {
  if (value < 0.5) { return "Low score: results may be incomplete or noisey" };
  if (value < 0.9) { return "Moderate score: results should be fairly reliable" };
  return "High score: results are likely very reliable";
}

export const WorkspaceItemCard: React.FC<WorkspaceItemCardProps> = ({
  item,
  selected,
  onToggleSelect,
  onView,
  onDelete,
  onRename,
}) => {
  const isCompound = item.kind === "compound"; // there are only two types: "compound" and "gene_cluster"

  // Tick every 15s so "X ago" updates
  const [, forceTick] = React.useState(0);
  React.useEffect(() => {
    const id = window.setInterval(() => {
      forceTick(n => n + 1);
    }, 5000);
    return () => { window.clearInterval(id); }
  }, [])

  const isQueued = item.status === "queued";
  const showSpinner = item.status === "processing";
  const isError = item.status === "error";
  const isDone = item.status === "done";

  // Rename state
  const [isEditing, setIsEditing] = React.useState(false);
  const [draftName, setDraftName] = React.useState(item.name);

  // Keep draft in sync when item.name changes from server
  React.useEffect(() => {
    if (!isEditing) {
      setDraftName(item.name);
    }
  }, [item.name, isEditing]);

  const startEditing = (e: React.MouseEvent) => {
    e.stopPropagation();
    setDraftName(item.name);
    setIsEditing(true);
  }

  const cancelEditing = () => {
    setDraftName(item.name);
    setIsEditing(false);
  }

  const commitEditing = () => {
    const trimmed = draftName.trim();
    // No-op if unchanged or empty
    if (!trimmed || trimmed === item.name) {
      setIsEditing(false);
      setDraftName(item.name);
      return;
    }
    onRename(item.id, trimmed);
    setIsEditing(false);
  }

  const handleNameKeyDown = (e: React.KeyboardEvent<HTMLInputElement>) => {
    if (e.key === "Enter") {
      e.preventDefault();
      commitEditing();
    } else if (e.key === "Escape") {
      e.preventDefault();
      cancelEditing();
    }
  }

  return (
    <Stack
      onClick={() => onToggleSelect(item.id)}
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
        sx={(theme) => {
          const t = theme.vars || theme;
          return {
            display: "flex",
            alignItems: "flex-start",
            justifyContent: "space-between",
            gap: 1.5,
          }
        }}
      >
        <Stack direction="row" spacing={1.5} alignItems="center">
          <Checkbox
            size="small"
            checked={selected}
            onClick={(e) => {
              e.stopPropagation();
              onToggleSelect(item.id);
            }}
          />

          {isCompound ? (
            <ScienceIcon fontSize="small" />
          ) : (
            <BiotechIcon fontSize="small" />
          )}

          <Stack direction="column" spacing={0.5}>
            <Stack direction="row" spacing={0.5} alignItems="center" width="400px">
              <Stack direction="row" alignItems="center" spacing={1}>
                <ItemKindChip itemKind={item.kind} />
                {isEditing ? (
                  <TextField
                    size="small"
                    variant="standard"
                    value={draftName}
                    autoFocus
                    onClick={(e) => e.stopPropagation()}
                    onChange={(e) => setDraftName(e.target.value)}
                    onKeyDown={handleNameKeyDown}
                    onBlur={commitEditing}
                    inputProps={{
                      style: { fontSize: "0.875rem", fontWeight: 500 },
                    }}
                  />
                ) : (
                  <>
                    <Typography
                      variant="body2"
                      fontWeight={500}
                      noWrap
                      sx={{ maxWidth: 220 }}
                    >
                      {item.name}
                    </Typography>
                    {isDone && (
                      <EditIcon
                        onClick={startEditing}
                        sx={{
                          fontSize: 16,
                          cursor: "pointer",
                          ml: 0.5
                        }}
                      />
                    )}
                  </>
                )}
              </Stack>
            </Stack>

            <Typography variant="caption" color="text.secondary">
              Status updated {formatUpdatedAgo(item.updatedAt)}
            </Typography>
          </Stack>
        </Stack>

        <Stack direction="row" spacing={1} alignItems="center">
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
            onClick={(e) => {
              e.stopPropagation();
              onView(item.id);
            }}
            disabled={!isDone}
          >
            <VisibilityIcon fontSize="small" />
          </IconButton>

          <IconButton
            size="small"
            onClick={(e) => {
              e.stopPropagation();
              onDelete(item.id);
            }}
          >
            <DeleteIcon fontSize="small" />
          </IconButton>
        </Stack>
      </Box>

      <Box>
        {isDone && (
          <Stack spacing={0.5} sx={{ pl: 8 }}>
            {item.retrofingerprints.map((fp, idx) => (
              <Stack
                key={fp.id}
                direction="row"
                spacing={1}
                alignItems="center"
              >
                <Typography
                  variant="caption"
                  color="text.secondary"
                  sx={{
                    lineHeight: 1,
                    minWidth: 60,
                    textAlign: "left",
                  }}
                >
                  {`Readout ${idx + 1}`}
                </Typography>
                <ScoreBar
                  key={fp.id}
                  score={fp.score}
                  getTooltipTitle={getScoreTooltip}
                  getScoreColor={getScoreColor}
                  width={200}
                  height={10}
                  tooltipPosition="right"
                />
              </Stack>
            ))}
          </Stack>
        )}
      </Box>
    </Stack>
  )
}