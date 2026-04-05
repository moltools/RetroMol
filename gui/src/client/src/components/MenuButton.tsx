import React from "react";
import Badge, { badgeClasses } from "@mui/material/Badge";
import IconButton, { IconButtonProps } from "@mui/material/IconButton";
import Tooltip from "@mui/material/Tooltip";

export interface MenuButtonProps extends IconButtonProps {
  showBadge?: boolean;
  tooltip?: string;
}

export const MenuButton: React.FC<MenuButtonProps> = ({ showBadge = false, tooltip = "", ...props }) => {
  return (
    <Badge
      color="error"
      variant="dot"
      invisible={!showBadge}
      sx={{ [`& .${badgeClasses.badge}`]: { right: 2, top: 2 } }}
    >
      <Tooltip title={tooltip} placement="bottom" arrow>
      <IconButton size="small" {...props} />
      </Tooltip>
    </Badge>
  )
}