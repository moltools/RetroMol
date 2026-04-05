import React from "react";
import Stack from "@mui/material/Stack";
import { NavbarBreadcrumbs } from "../NavbarBreadcrumbs";
import { NotificationDrawer } from "../NotificationDrawer";
import { WorkspaceControls } from "./WorkspaceControls";

export const WorkspaceHeader: React.FC = () => {
  const [drawerOpen, setDrawerOpen] = React.useState(false);

  // Function to handle the opening of the notification drawer
  const handleDrawerOpen = () => {
    setDrawerOpen(true);
  }

  // Function to handle the closing of notification drawer
  const handleDrawerClose = () => {
    setDrawerOpen(false);
  }

  return (
    <>
      <Stack
        direction="row"
        sx={{
          display: { xs: "none", md: "flex" },
          width: "100%",
          alignItems: { xs: "flex-start", md: "center" },
          justifyContent: "space-between",
          maxWidth: { sm: "100%", md: "1700px" },
          pt: 1.5,
        }}
        spacing={2}
      >
        <NavbarBreadcrumbs />
        <WorkspaceControls handleDrawerOpen={handleDrawerOpen} />
      </Stack>
      <NotificationDrawer open={drawerOpen} handleClose={handleDrawerClose} />
    </>
  )
}