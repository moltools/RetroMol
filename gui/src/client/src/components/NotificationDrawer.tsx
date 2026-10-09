import React, { useEffect } from "react";
import Box from "@mui/material/Box";
import MuiDrawer, { drawerClasses } from "@mui/material/Drawer";
import ListItem from "@mui/material/ListItem";
import ListItemText from "@mui/material/ListItemText";
import Stack from "@mui/material/Stack";
import Typography from "@mui/material/Typography";
import CloseIcon from "@mui/icons-material/Close";
import ClearAllIcon from "@mui/icons-material/ClearAll";
import { styled, useTheme } from "@mui/material/styles";
import { useNotifications } from "./NotificationProvider";

// Custom styling for the drawer paper (the visible panel). Deliberately does NOT set
// width/etc. on the Drawer's own root -- that root is the fixed, full-viewport-covering
// container MUI uses internally to compute how far off-screen the panel starts before
// sliding in; shrinking it to the panel's own width corrupts that calculation and makes
// the panel appear to enter from the left instead of the right.
const drawerWidth = 400;
const Drawer = styled(MuiDrawer)({
  [`& .${drawerClasses.paper}`]: {
    width: drawerWidth,
    boxSizing: "border-box",
  },
})

interface NotificationDrawerProps {
  open: boolean;
  handleClose: () => void;
}

export const NotificationDrawer: React.FC<NotificationDrawerProps> = ({ open, handleClose }) => {
  const { notifications, markAllAsRead, clearNotifications } = useNotifications();
  const theme = useTheme();

  // When the drawer opens, mark all notifications as read
  useEffect(() => {
    if (open) { markAllAsRead(); }
  }, [open, markAllAsRead]);

  return (
    <Drawer
    anchor="right"
    open={open}
    onClose={handleClose}
  >
      <Box
        sx={{
          width: drawerWidth,
          height: "100%",
          bgcolor: "background.paper",
          borderLeft: "1px solid",
          borderColor: "divider",
          pt: 1.5,
          overflowY: "auto",
          display: "flex",
        }}
      >
        <Stack spacing={2} width="100%">
          <Box
            sx={{
              display: "flex",
              alignItems: "center",
              justifyContent: "space-between",
              px: 3,
              pb: 0,
            }}
          >
            <Typography variant="subtitle1">Notifications</Typography>
            <Box>
              <ClearAllIcon
                onClick={clearNotifications}
                aria-label="Clear all notifications"
                sx={{
                  mr: 3,
                  cursor: "pointer",
                }}
              />
              <CloseIcon
                onClick={handleClose}
                aria-label="Close notifications"
                sx={{
                  position: "absolute",
                  right: 16,
                  cursor: "pointer",
                }}
              />
            </Box>
          </Box>

          {notifications.length === 0 ? (
            <Typography
              variant="body2"
              sx={{
                px: 3,
                pt: "4px",
                pb: 2,
              }}
            >
              No notifications to show.
            </Typography>
          ) : (
            <Stack
              sx={{
                overflowY: "auto",
                px: 3,
                pt: "4px",
                pb: 2,
              }}
            >
              {notifications.slice(0).reverse().map((notification) => (
                <ListItem
                  key={crypto.randomUUID()}
                  divider
                  sx={{
                    borderLeft: `4px solid ${theme.palette[notification.level].main}`,
                    backgroundColor: "transparent",
                    transition: "background-color 0.3s ease",
                  }}
                >
                  <ListItemText
                    primary={ new Date(notification.timestamp).toLocaleString() }
                    secondary={
                      <Typography
                        variant="subtitle2"
                        sx={{ fontWeight: notification.isNew ? 600 : 400 }}
                      >
                        {notification.content}
                      </Typography>
                    }
                    secondaryTypographyProps={{ component: "div" }}  // render wrapper as div instead of p
                  />
                </ListItem>
              ))}
            </Stack>
          )}
        </Stack>
      </Box>
    </Drawer>
  )
}