import React from "react";
import { Box, IconButton, IconButtonOwnProps, Menu, MenuItem, Tooltip } from "@mui/material";
import { DarkMode as DarkModeIcon, LightMode as LightModeIcon } from "@mui/icons-material";
import { useColorScheme } from "@mui/material/styles";

export default function ColorModeIconDropdown(props: IconButtonOwnProps) {
  const { mode, systemMode, setMode } = useColorScheme();
  const [anchorEl, setAnchorEl] = React.useState<null | HTMLElement>(null);
  const open = Boolean(anchorEl);

  // Handle the opening and closing of the menu
  const handleClick = (event: React.MouseEvent<HTMLElement>) => {
    setAnchorEl(event.currentTarget);
  }

  // Handle the closing of the menu
  const handleClose = () => {
    setAnchorEl(null);
  }

  // Handle the selection of a mode
  const handleMode = (targetMode: "system" | "light" | "dark") => () => {
    setMode(targetMode);
    handleClose();
  }

  if (!mode) {
    return (
      <Box
        data-screenshot="toggle-mode"
        sx={(theme) => ({
          verticalAlign: "bottom",
          display: "inline-flex",
          width: "2.25rem",
          height: "2.25rem",
          borderRadius: (theme.vars || theme).shape.borderRadius,
          border: "1px solid",
          borderColor: (theme.vars || theme).palette.divider,
        })}
      />
    );
  }

  const resolvedMode = (systemMode || mode) as "light" | "dark";
  const icon = { light: <LightModeIcon />, dark: <DarkModeIcon /> }[resolvedMode];

  return (
    <React.Fragment>
      <Tooltip title={"Open color mode menu"} placement="bottom" arrow>
        <IconButton
          data-screenshot="toggle-mode"
          onClick={handleClick}
          disableRipple
          size="small"
          aria-controls={open ? "color-scheme-menu" : undefined}
          aria-haspopup="true"
          aria-expanded={open ? "true" : undefined}
          {...props}
        >
          {icon}
        </IconButton>
      </Tooltip>
      <Menu
        anchorEl={anchorEl}
        id="account-menu"
        open={open}
        onClose={handleClose}
        onClick={handleClose}
        slotProps={{
          paper: {
            variant: "outlined",
            elevation: 0,
            sx: { my: "4px", p: "8px", width: "max-content", minWidth: "200px" },
          },
        }}
        transformOrigin={{ horizontal: "right", vertical: "top" }}
        anchorOrigin={{ horizontal: "right", vertical: "bottom" }}
      >
        <MenuItem selected={mode === "system"} onClick={handleMode("system")} sx={{ mb: 1 }}>
          System
        </MenuItem>
        <MenuItem selected={mode === "light"} onClick={handleMode("light")} sx={{ mb: 1 }}>
          Light
        </MenuItem>
        <MenuItem selected={mode === "dark"} onClick={handleMode("dark")}>
          Dark
        </MenuItem>
      </Menu>
    </React.Fragment>
  );
}
