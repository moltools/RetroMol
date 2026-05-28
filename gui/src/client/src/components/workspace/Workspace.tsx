import React from "react";
import Box from "@mui/material/Box";
import Stack from "@mui/material/Stack";
import Fade from "@mui/material/Fade";
import { alpha } from "@mui/material/styles";
import { Routes, Route, useNavigate } from "react-router-dom";
import { useNotifications } from "../NotificationProvider";
import { useOverlay } from "../OverlayProvider";
import { Session } from "../../features/session/types";
import { getSession, refreshSession } from "../../features/session/api";
import { WorkspaceNavbar } from "./WorkspaceNavbar";
import { WorkspaceSideMenu } from "./WorkspaceSideMenu";
import { WorkspaceHeader } from "./WorkspaceHeader";
import { WorkspaceHome } from "./WorkspaceHome";
import { WorkspaceUpload } from "./WorkspaceUpload";
import { WorkspaceDiscovery } from "./WorkspaceDiscovery";
// import { WorkspaceEnrichment } from "./tabs/enrichment/WorkspaceEnrichment";

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

  // SSE: refresh session when server says something changed
  React.useEffect(() => {
    if (!session?.sessionId) return;

    let alive = true;
    let refreshTimer: number | null = null;

    const scheduleRefresh = () => {
      if (!alive) return;
      if (refreshTimer !== null) return;

      refreshTimer = window.setTimeout(() => {
        refreshTimer = null;

        refreshSession(session.sessionId)
          .then((fresh) => {
            if (!alive) return;
            setSession(fresh);
          })
          .catch((err) => {
            const msg = err instanceof Error ? err.message : String(err);
            pushNotification(`Failed to refresh session: ${msg}`, "error");
          });
      }, 250);
    };

    const SSE_BASE = process.env.REACT_APP_SSE_BASE ?? "";
    const es = new EventSource(
      `${SSE_BASE}/api/sessionEvents?sessionId=${encodeURIComponent(session.sessionId)}`
    );

    // Attach scheduleRefresh to all relevant events
    es.addEventListener("hello", scheduleRefresh);
    es.addEventListener("item_updated", scheduleRefresh);
    es.addEventListener("session_merged", scheduleRefresh);

    es.onopen = () => {
      // EventSource retries automatically; avoid spamming notifications
    };

    return () => {
      alive = false;
      if (refreshTimer !== null) window.clearTimeout(refreshTimer);
      es.close();
    };
  }, [session?.sessionId, pushNotification]);

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
                  {/*<Route path="enrichment" element={<WorkspaceEnrichment session={session} setSession={setSession} />} />*/}
                    <Route path="enrichment" element={<div>Analysis currently available. Check back later.</div>} />
                </Routes>
              </Box>
            </Fade>
          )}

        </Stack>
      </Box>
    </Box>
  )
};