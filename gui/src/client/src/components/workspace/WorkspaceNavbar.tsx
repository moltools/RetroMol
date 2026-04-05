import React from "react";
import { styled } from "@mui/material/styles";
import AppBar from "@mui/material/AppBar";
import MuiToolbar from "@mui/material/Toolbar";
import Stack from "@mui/material/Stack";
import Typography from "@mui/material/Typography";
import { tabsClasses } from "@mui/material/Tabs";
import MenuRoundedIcon from "@mui/icons-material/MenuRounded";
import { WorkspaceSideMenuMobile } from "./WorkspaceSideMenuMobile";
import { MenuButton } from "../MenuButton";

const Toolbar = styled(MuiToolbar)({
  width: "100%",
  padding: "12px",
  display: "flex",
  flexDirection: "column",
  alignItems: "start",
  justifyContent: "center",
  gap: "12px",
  flexShrink: 0,
  [`& ${tabsClasses.flexContainer}`]: { gap: "8px", p: "8px", pb: 0 },
});

export const WorkspaceNavbar: React.FC = () => {
  const [open, setOpen] = React.useState(false);

  // Toggle the mobile menu drawer
  const toggleDrawer = (newOpen: boolean) => () => {
    setOpen(newOpen);
  }

  return (
    <AppBar
      position="fixed"
      sx={{
        display: { xs: "auto", md: "none" },
        boxShadow: 0,
        bgcolor: "background.paper",
        backgroundImage: "none",
        borderBottom: "1px solid",
        borderColor: "divider",
        top: "var(--template-frame-height, 0px)",
      }}
    >
      <Toolbar variant="regular">
        <Stack direction="row" sx={{ alignItems: "center", flexGrow: 1, width: "100%", gap: 1 }}>
          <Stack direction="row" spacing={1} sx={{ justifyContent: "center", mr: "auto" }}>
            <Typography variant="h4" component="h1" sx={{ color: "text.primary" }}>
            </Typography>
          </Stack>
          <MenuButton aria-label="menu" onClick={toggleDrawer(true)}>
            <MenuRoundedIcon />
          </MenuButton>
          <WorkspaceSideMenuMobile open={open} toggleDrawer={toggleDrawer} />
        </Stack>
      </Toolbar>
    </AppBar>
  )
}