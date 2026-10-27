import React from "react";
import Box from "@mui/material/Box";
import Stack from "@mui/material/Stack";
import Fade from "@mui/material/Fade";
import { alpha } from "@mui/material/styles";
import { Routes, Route, useNavigate } from "react-router-dom";
import { useNotifications } from "../NotificationProvider";
import { useOverlay } from "../OverlayProvider";
import { Session } from "../../features/session/types";
import { getSession } from "../../features/session/api";
import { useSessionPolling } from "../../features/session/useSessionPolling";
import { WorkspaceNavbar } from "./WorkspaceNavbar";
import { WorkspaceSideMenu } from "./WorkspaceSideMenu";
import { WorkspaceHeader } from "./WorkspaceHeader";
import { WorkspaceHome } from "./WorkspaceHome";
import { WorkspaceUpload } from "./WorkspaceUpload";
import { WorkspaceDiscovery } from "./WorkspaceDiscovery";
import { WorkspaceRules } from "./WorkspaceRules";
import { WorkspaceGenerate } from "./WorkspaceGenerate";

export const Workspace: React.FC = () => {
  const { showOverlay, hideOverlay } = useOverlay();
  const { pushNotification } = useNotifications();
  const navigate = useNavigate();

  const [loading, setLoading] = React.useState<boolean>(true);
  const [session, setSession] = React.useState<Session | null>(null);

  // Load session on mount
  React.useEffect(() => {
    let alive = true;
    setLoading(true);

    getSession()
      .then((sess) => {
        if (!alive) return;
        setSession(sess);
      })
      .catch((err) => {
        console.error("Error loading session:", err);
        navigate("/notfound");
      })
      .finally(() => {
        if (!alive) return;
        setLoading(false);
      });

    return () => { alive = false; };
  }, [navigate]);

  // Overlay follows loading state
  React.useEffect(() => {
    if (loading) showOverlay();
    else hideOverlay();
  }, [loading, showOverlay, hideOverlay]);

  const reportRefreshError = React.useCallback(
    (message: string) => pushNotification(message, "error"),
    [pushNotification],
  );
  useSessionPolling(session, setSession, reportRefreshError);

  // Determine what to show
  const showContent = !!session && !loading;

  if (!session && !loading) return null;

  return (
    <Box sx={{ display: "flex"}}>
      <WorkspaceNavbar />
      <WorkspaceSideMenu />
      <Box
        component="main"
        sx={theme => ({
          flexGrow: 1,
          backgroundColor: theme.vars
            ? `rgba(${theme.vars.palette.background.defaultChannel} /  1)`
            : alpha(theme.palette.background.default, 1),
          overflow: "auto",
        })}
      >
        <Stack
          spacing={2}
          sx={{ alignItems: "center", mx: 3, pb: 5, mt: { xs: 8, md: 0 } }}
        >
          <WorkspaceHeader />

          {/* Actual workspace content fades in once session is ready */}
          {session && (
            <Fade in={showContent} timeout={200} unmountOnExit>
              <Box sx={{ width: "100%" }}>
                <Routes>
                  <Route index element={<WorkspaceHome />} />
                  <Route path="upload" element={<WorkspaceUpload session={session} setSession={setSession} />} />
                  <Route path="discovery" element={<WorkspaceDiscovery session={session} setSession={setSession} />} />
                  <Route path="generate" element={<WorkspaceGenerate />} />
                  <Route path="rules" element={<WorkspaceRules />} />
                </Routes>
              </Box>
            </Fade>
          )}

        </Stack>
      </Box>
    </Box>
  )
};