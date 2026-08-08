import React from "react";
import CircularProgress from "@mui/material/CircularProgress";
import IconButton from "@mui/material/IconButton";
import Tooltip from "@mui/material/Tooltip";
import DownloadIcon from "@mui/icons-material/Download";
import { useNotifications } from "./NotificationProvider";
import { exportElementAsPng } from "./exportImage";

// Small "export this view as a PNG" affordance. Captures whatever is currently
// rendered inside `targetRef`, including any live highlight state, exactly as
// shown on screen.
export function ExportImageButton({
  targetRef,
  filename,
  label = "Download this view as a PNG",
}: {
  targetRef: React.RefObject<HTMLElement>;
  filename: string;
  label?: string;
}) {
  const { pushNotification } = useNotifications();
  const [exporting, setExporting] = React.useState(false);

  const handleExport = async () => {
    const node = targetRef.current;
    if (!node) return;

    setExporting(true);
    try {
      await exportElementAsPng(node, filename);
    } catch (err) {
      const msg = err instanceof Error ? err.message : String(err);
      pushNotification(`Failed to export image: ${msg}`, "error");
    } finally {
      setExporting(false);
    }
  };

  return (
    <Tooltip title={label} arrow>
      <span>
        <IconButton size="small" disabled={exporting} onClick={handleExport}>
          {exporting ? <CircularProgress size={18} /> : <DownloadIcon fontSize="small" />}
        </IconButton>
      </span>
    </Tooltip>
  );
}
