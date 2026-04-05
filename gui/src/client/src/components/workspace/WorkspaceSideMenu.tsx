import React from "react";
import Box from "@mui/material/Box";
import MuiDrawer, { drawerClasses } from "@mui/material/Drawer";
import { styled } from "@mui/material/styles";
import { MenuContent } from "../MenuContent";
import { BuildVersion } from "../BuildVersion";

const drawerWidth = 240;

const Drawer = styled(MuiDrawer)({
  width: drawerWidth,
  flexShrink: 0,
  boxSizing: "border-box",
  mt: 10,
  [`& .${drawerClasses.paper}`]: {
    width: drawerWidth,
    boxSizing: "border-box",
  },
})

export const WorkspaceSideMenu: React.FC = () => {
  return (
    <Drawer
      variant="permanent"
      sx={{
        display: { xs: "none", md: "block" },
        [`& .${drawerClasses.paper}`]: {
          backgroundColor: "background.paper",
        },
      }}
    >
      <Box
        sx={{
          overflow: "auto",
          height: "100%",
          display: "flex",
          flexDirection: "column",
        }}
      >
        <MenuContent />
        <BuildVersion />
      </Box>
    </Drawer>
  )
}