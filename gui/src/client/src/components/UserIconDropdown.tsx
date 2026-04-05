import React from "react";
import { useNavigate } from "react-router-dom";
import Button from "@mui/material/Button";
import Dialog from "@mui/material/Dialog";
import DialogActions from "@mui/material/DialogActions";
import DialogContent from "@mui/material/DialogContent";
import DialogContentText from "@mui/material/DialogContentText";
import DialogTitle from "@mui/material/DialogTitle";
import IconButton, { IconButtonOwnProps } from "@mui/material/IconButton";
import Menu from "@mui/material/Menu";
import MenuItem from "@mui/material/MenuItem";
import Tooltip from "@mui/material/Tooltip";
import AccountCircleIcon from "@mui/icons-material/AccountCircle";
import { useNotifications } from "./NotificationProvider";
import { deleteSession } from "../features/session/api";
import { deleteCookie, getCookie } from "../features/session/utils";

export const UserIconDropdown: React.FC<IconButtonOwnProps> = (props) => {
  const { pushNotification } = useNotifications();
  const navigate = useNavigate();

  const [anchorEl, setAnchorEl] = React.useState<null | HTMLElement>(null);
  const open = Boolean(anchorEl);

  const [deleteDialogOpen, setDeleteDialogOpen] = React.useState(false);
  const [logoutDialogOpen, setLogoutDialogOpen] = React.useState(false);

  const handleCopy = async () => {
    // Get sessionId from cookies
    const sessionId = getCookie("sessionId");
    if (!sessionId) {
      pushNotification("Session ID not found in cookies", "error");
      return;
    };

    try {
      // Copy sessionId to clipboard
      await navigator.clipboard.writeText(sessionId);
      pushNotification("Session ID copied to clipboard", "success");
    } catch (err) {
      pushNotification("Failed to copy session ID", "error");
    };
  };

  // Remove session ID from cookies and redirect to home page
  const confirmLogout = () => {
    deleteCookie("sessionId");
    setLogoutDialogOpen(false);
    navigate("/");
  }

  // Call API to delete session after confirmation
  const confirmDeleteSession = async () => {
    try {
      await deleteSession();
      deleteCookie("sessionId");
      setDeleteDialogOpen(false);
      navigate("/");
    } catch (err) {
      console.error("Error deleting session:", err);
      const msg = err instanceof Error ? err.message : String(err);
      pushNotification(`Failed to delete session: ${msg}`, "error");
    }
  }

  return (
    <>
    <React.Fragment>
      <Tooltip title="Open user menu" placement="bottom" arrow>
        <IconButton
          onClick={(event => { setAnchorEl(event.currentTarget) })}
          disableRipple
          size="small"
          {...props}
        >
          <AccountCircleIcon />
        </IconButton>
      </Tooltip>
      <Menu
        anchorEl={anchorEl}
        id="account-menu"
        open={open}
        onClose={() => { setAnchorEl(null) }}
        onClick={() => { setAnchorEl(null) }}
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
        <MenuItem onClick={handleCopy} sx={{ mb: 1 }}>
          Copy session ID
        </MenuItem>
        <MenuItem onClick={() => { setLogoutDialogOpen(true) }} sx={{ mb: 1 }}>
          Logout session
        </MenuItem>
        <MenuItem onClick={() => { setDeleteDialogOpen(true) }}>
          Delete session
        </MenuItem>
      </Menu>
    </React.Fragment>

    {/* Delete confirmation dialog */}
    <Dialog
      open={deleteDialogOpen}
      onClose={() => setDeleteDialogOpen(false)}
      container={() => document.body}
      fullWidth
      maxWidth="sm"
    >
      <DialogTitle>Delete session</DialogTitle>
      <DialogContent>
        <DialogContentText>
          Do you really want to delete your session?<br /><br />
          All sessions will be automatically deleted <b>seven days</b> after session inactivity.
          Deleting your session will remove all your uploaded data immediately. You will not be able to retrieve it later.
        </DialogContentText>
      </DialogContent>
      <DialogActions>
        <Button onClick={() => setDeleteDialogOpen(false)} color="primary">
          Cancel
        </Button>
        <Button onClick={confirmDeleteSession} color="error">
          Delete
        </Button>
      </DialogActions>
    </Dialog>

    {/* Logout confirmation dialog */}
    <Dialog
      open={logoutDialogOpen}
      onClose={() => setLogoutDialogOpen(false)}
      container={() => document.body}
      fullWidth
      maxWidth="sm"
    >
      <DialogTitle>Logout session</DialogTitle>
      <DialogContent>
        <DialogContentText>
          Do you really want to log out of your session?<br /><br />
          Make sure to copy your session ID before logging out. You will not be able to retrieve your session later if you lose it. Administrators will not be able to access your session or retrieve your session ID for you.
        </DialogContentText>
      </DialogContent>
      <DialogActions>
        <Button onClick={() => setLogoutDialogOpen(false)} color="primary">
          Cancel
        </Button>
        <Button onClick={confirmLogout} color="primary">
          Logout
        </Button>
      </DialogActions>
    </Dialog>
    </>
  )
}