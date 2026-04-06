import React from "react";
import Box from "@mui/material/Box";
import Stack from "@mui/material/Stack";
import Tooltip from "@mui/material/Tooltip";
import Typography from "@mui/material/Typography";
import { alpha } from "@mui/material/styles";

type ScoreBarProps = {
  score: number; // value between 0 and 1
  getTooltipTitle: (value: number) => string;
  getScoreColor: (theme: any, value: number) => string;
  tooltipPosition?: "top" | "bottom" | "left" | "right";
  width?: number | string;
  height?: number | string;
}

export const ScoreBar: React.FC<ScoreBarProps> = ({
  score,
  getTooltipTitle,
  getScoreColor,
  tooltipPosition = "top",
  width = "100%",
  height = 16
}) => {
  const value = Math.max(0, Math.min(1, score));

  return (
    <Box sx={{ mt: 0.5 }}>
      <Tooltip title={getTooltipTitle(value)} placement={tooltipPosition} arrow>
        <Stack direction="row" spacing={0.5} alignItems="center">
          <Box
            sx={(theme) => {
              const t = theme.vars || theme;
              return {
                height: height,
                width: width,
                flexShrink: 0,
                borderRadius: 999,
                overflow: "hidden",
                backgroundColor: alpha("#000000", 0.1),
                ...theme.applyStyles("dark", { backgroundColor: alpha("#ffffff", 0.1) })
            }}}
          >
            <Box
              sx={(theme) => {
                const t = theme.vars || theme;
                const barColor = getScoreColor(t, value);
                return {
                  width: `${value * 100}%`,
                  height: "100%",
                  backgroundColor: barColor,
                  transition: "width 200ms ease-out",
                };
              }}
            />
          </Box>
          <Typography variant="caption" sx={{ minWidth: 30, textAlign: "right" }}>
            {(value * 100).toFixed(1)}%
          </Typography>
        </Stack>
      </Tooltip>
    </Box>
  )
}