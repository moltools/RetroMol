import React from "react";
import Tooltip from "@mui/material/Tooltip";
import ContentCopyIcon from "@mui/icons-material/ContentCopy";
import { useNotifications } from "./NotificationProvider";
import { MinimalIconButton } from "./MinimalIconButton";

// Copies `text` to the clipboard on click, with a toast for success/failure --
// same pattern UserIconDropdown's session-id copy uses, pulled out for reuse
// wherever a raw string (a SMILES, a reaction SMARTS, ...) needs a copy affordance.
export function CopyIconButton({ text, label = "value" }: { text: string; label?: string }) {
  const { pushNotification } = useNotifications();

  const handleCopy = async (event: React.MouseEvent) => {
    event.stopPropagation();
    try {
      await navigator.clipboard.writeText(text);
      pushNotification(`Copied ${label} to clipboard`, "success");
    } catch (err) {
      pushNotification(`Failed to copy ${label}`, "error");
    }
  };

  return (
    <Tooltip title={`Copy ${label}`} arrow>
      <MinimalIconButton onClick={handleCopy}>
        <ContentCopyIcon fontSize="inherit" />
      </MinimalIconButton>
    </Tooltip>
  );
}
