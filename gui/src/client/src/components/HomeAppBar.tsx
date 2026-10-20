import React from "react";
import { styled, alpha } from "@mui/material/styles";
import AppBar from "@mui/material/AppBar";
import Box from "@mui/material/Box";
import Button from "@mui/material/Button";
import Container from "@mui/material/Container";
import Divider from "@mui/material/Divider";
import Drawer from "@mui/material/Drawer";
import IconButton from "@mui/material/IconButton";
import InputAdornment from "@mui/material/InputAdornment";
import InputLabel from "@mui/material/InputLabel";
import MenuItem from "@mui/material/MenuItem";
import TextField from "@mui/material/TextField";
import Toolbar from "@mui/material/Toolbar";
import MenuIcon from "@mui/icons-material/Menu";
import CloseRoundedIcon from "@mui/icons-material/CloseRounded";
import visuallyHidden from "@mui/utils/visuallyHidden";
import ColorModeIconDropdown from "../theme/ColorModeIconDropdown";
import { useNavigate } from "react-router-dom";
import { createSession, getSession } from "../features/session/api";
import { createCookie } from "../features/session/utils";

// Custom styling for the toolbar
const StyledToolbar = styled(Toolbar)(({ theme }) => ({
  display: "flex",
  alignItems: "center",
  justifyContent: "space-between",
  flexShrink: 0,
  borderRadius: `calc(${theme.shape.borderRadius}px + 8px)`,
  backdropFilter: "blur(24px)",
  border: "1px solid",
  borderColor: (theme.vars || theme).palette.divider,
  backgroundColor: theme.vars
    ? `rgba(${theme.vars.palette.background.defaultChannel} / 0.4)`
    : alpha(theme.palette.background.default, 0.4),
  boxShadow: (theme.vars || theme).shadows[1],
  padding: "8px 12px",
}))


// Component for session ID input
function RetrieveSession() {
  const navigate = useNavigate();
  const [sessionId, setSessionId] = React.useState<string>("");

  // Function that checks if the session ID is valid, and if so, redirects to the dashboard
  const handleRetrieveSession = async () => {
    // Validate: empty sessionId is immediately invalid
    if (!sessionId.trim()) {
      navigate("/notfound");
      return;
    }

    // Try to retrieve the session from the backend
    try {
      await getSession(sessionId);
      createCookie("sessionId", sessionId);
      navigate("/dashboard");
    } catch (err) {
      console.error("error retrieving session:", err);
      navigate("/notfound");
    }
  }

  return (
    <>
      <InputLabel htmlFor="session-id-input" sx={visuallyHidden}>
        Session ID
      </InputLabel>
      <TextField
        id="session-id-input"
        hiddenLabel
        size="small"
        variant="outlined"
        aria-label="Enter your session ID"
        placeholder="Your session ID"
        fullWidth
        value={sessionId}
        onChange={(e) => setSessionId(e.target.value)}
        slotProps={{
          input: {
            autoComplete: "off",
            "aria-label": "Enter your session ID",
            endAdornment: (
              <InputAdornment position="end">
                <Button
                  color="primary"
                  variant="contained"
                  size="small"
                  onClick={handleRetrieveSession}
                  sx={{ minWidth: "fit-content", height: "24px" }}
                >
                  Retrieve
                </Button>
              </InputAdornment>
            )
          },
        }}
      />
    </>
  )
}

export default function HomeAppBar() {
  const navigate = useNavigate();
  const [open, setOpen] = React.useState(false);

  // Function that opens FAQ page in a new tab
  const handleOpenFAQ = () => {
    window.open(
      "https://github.com/moltools/RetroMol/issues",
      "_blank",
      "noopener,noreferrer"
    );
  }

  // Function to handle the opening and closing of the drawer
  const toggleDrawer = (newOpen: boolean) => () => {
    setOpen(newOpen);
  }

  // Function to handle the creation of a new session
  const handleCreateSession = async (event: React.FormEvent) => {
    event.preventDefault();

    try {
      const sessionId = await createSession();
      createCookie("sessionId", sessionId);
      navigate("/dashboard");
    } catch (err) {
      const status = (err as any)?.status;
      if (status === 503) {
        navigate("/maxsessions");
      } else {
        navigate("/oops");
      }
    }
  }

  return (
    <AppBar
      position="fixed"
      enableColorOnDark
      sx={{
        boxShadow: 0,
        bgcolor: "transparent",
        backgroundImage: "none",
        mt: "calc(var(--template-frame-height, 0px) + 28px)",
      }}
    >
      <Container maxWidth="lg">
        <StyledToolbar variant="dense" disableGutters>
          <Box sx={{ flexGrow: 1, display: "flex", alignItems: "center", px: 0 }}>
            <Box sx={{ display: { xs: "none", md: "flex" } }}>
              <Button
                variant="text"
                color="info"
                size="small"
                sx={{ minWidth: 0 }}
                onClick={handleOpenFAQ}
              >
                FAQ
              </Button>
            </Box>
          </Box>
          <Box
            sx={{
              display: { xs: "none", md: "flex" },
              gap: 1,
              alignItems: "center",
            }}
          >
            <RetrieveSession />
            <Button
              color="primary"
              variant="contained"
              size="small"
              onClick={handleCreateSession}
              sx={{ minWidth: "fit-content" }}
            >
              New session
            </Button>
            <ColorModeIconDropdown />
          </Box>
          <Box sx={{ display: { xs: "flex", md: "none" }, gap: 1 }}>
            <ColorModeIconDropdown size="medium" />
            <IconButton aria-label="Menu button" onClick={toggleDrawer(true)}>
              <MenuIcon />
            </IconButton>
            <Drawer anchor="top" open={open} onClose={toggleDrawer(false)}>
              <Box sx={{ p: 2, backgroundColor: "background.default" }}>
                <Box sx={{ display: "flex", justifyContent: "flex-end" }}>
                  <IconButton onClick={toggleDrawer(false)}>
                    <CloseRoundedIcon />
                  </IconButton>
                </Box>

                <MenuItem onClick={() => { handleOpenFAQ(); setOpen(false); }}>
                  FAQ
                </MenuItem>

                <Divider sx={{ my: 3 }} />

                <MenuItem>
                  <Button
                    color="primary"
                    variant="contained"
                    fullWidth
                    onClick={handleCreateSession}
                  >
                    New session
                  </Button>
                </MenuItem>
                <MenuItem>
                  <RetrieveSession />
                </MenuItem>
              </Box>
            </Drawer>
          </Box>
        </StyledToolbar>
      </Container>
    </AppBar>
  )
}