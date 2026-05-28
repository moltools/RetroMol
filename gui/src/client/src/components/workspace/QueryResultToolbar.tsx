import React from "react";
import Box from "@mui/material/Box";
import Tooltip from "@mui/material/Tooltip";
import ZoomOutIcon from "@mui/icons-material/ZoomOut";
import ZoomInIcon from "@mui/icons-material/ZoomIn";
import RefreshIcon from "@mui/icons-material/Refresh";
import DownloadIcon from "@mui/icons-material/Download";
import ExchangeIcon from "@mui/icons-material/SwapHoriz";

type Props = {
  onZoomOut: () => void;
  onZoomIn: () => void;
  onZoomReset: () => void;
  onInvertMotifs: () => void;
  onDownloadSvg: () => void;
};

export const QueryResultToolbar: React.FC<Props> = ({
  onZoomOut,
  onZoomIn,
  onZoomReset,
  onInvertMotifs,
  onDownloadSvg,
}) => {
  return (
    <Box
      sx={{
        display: "flex",
        flexDirection: "row",
        justifyContent: "right",
        alignItems: "center",
        px: 2,
        py: 1,
      }}
    >
      <Box sx={{ display: "flex", alignItems: "center", gap: 1 }}>
        <Tooltip title="Zoom MSA out" arrow>
          <ZoomOutIcon onClick={onZoomOut} sx={{ cursor: "pointer" }} />
        </Tooltip>
        <Tooltip title="Zoom MSA in" arrow>
          <ZoomInIcon onClick={onZoomIn} sx={{ cursor: "pointer" }} />
        </Tooltip>
        <Tooltip title="Reset zoom" arrow>
          <RefreshIcon onClick={onZoomReset} sx={{ cursor: "pointer" }} />
        </Tooltip>
        <Tooltip title="Invert order of motifs" arrow>
          <ExchangeIcon onClick={onInvertMotifs} sx={{ cursor: "pointer" }} />
        </Tooltip>
        <Tooltip title="Download as SVG" arrow>
          <DownloadIcon onClick={onDownloadSvg} sx={{ cursor: "pointer" }} />
        </Tooltip>
      </Box>
    </Box>
  );
};