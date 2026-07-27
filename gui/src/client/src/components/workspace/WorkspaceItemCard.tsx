import React from "react";
import Box from "@mui/material/Box";
import Stack from "@mui/material/Stack";
import Typography from "@mui/material/Typography";
import Checkbox from "@mui/material/Checkbox";
import Chip from "@mui/material/Chip";
import IconButton from "@mui/material/IconButton";
import Tooltip from "@mui/material/Tooltip";
import DeleteIcon from "@mui/icons-material/Delete";
import ViewIcon from "@mui/icons-material/Visibility";
import CircularProgress from "@mui/material/CircularProgress";
import { Gauge } from "@mui/x-charts/Gauge";
import { SessionItem } from "../../features/session/types";
import { alpha } from "@mui/material/styles";
import type { Theme } from "@mui/material/styles";
import { DialogViewItem } from "./DialogViewItem";
import { useTick } from "../../hooks/useTick";

function getScoreColor(theme: Theme, value: number): string {
  const t = theme.vars || theme;
  if (value < 0.5) { return t.palette.error.main };
  if (value < 0.9) { return t.palette.warning.main };
  return t.palette.success.main;
};

type WorkspaceItemCardProps = {
  sessionId: string;
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
  sessionId,
  item,
  selected,
  disabled = false,
  onToggleSelect,
  onDelete,
}) => {
  const isCompound = item.kind === "compound"; // there are only two types: "compound" and "cluster"
  const itemScore = typeof item.score === "number" ? item.score : 0.0;

  const [openViewItem, setOpenViewItem] = React.useState(false);

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
              disabled={disabled || !isCompound}
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
          </Stack>
        </Box>
      </Stack>

      <DialogViewItem
        sessionId={sessionId}
        item={item}
        open={openViewItem}
        onClose={() => setOpenViewItem(false)}
      />
    </>
  );
};