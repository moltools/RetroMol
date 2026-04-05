import React from "react";
import Stack from "@mui/material/Stack";
import NotificationsRoundedIcon from "@mui/icons-material/NotificationsRounded";
import { MenuButton } from "../MenuButton";
import { useNotifications } from "../NotificationProvider";
import { UserIconDropdown } from "../UserIconDropdown";
import ColorModeIconDropdown from "../../theme/ColorModeIconDropdown";

interface WorksspaceControlsProps {
  handleDrawerOpen: () => void;
  sx?: object;
}

export const WorkspaceControls: React.FC<WorksspaceControlsProps> = ({ handleDrawerOpen, sx }) => {
  const { notifications } = useNotifications();

  // Check if there"s any new notification
  const hasNewNotifications = notifications.some((n) => n.isNew);

  return (
    <Stack direction="row" sx={{ gap: 1, ...sx }}>
      <MenuButton
        showBadge={hasNewNotifications}
        aria-label="Open notifications"
        onClick={handleDrawerOpen}
        tooltip="Open notifications drawer"
      >
          <NotificationsRoundedIcon />
      </MenuButton>
      <ColorModeIconDropdown />
      <UserIconDropdown />
    </Stack>
  )
}