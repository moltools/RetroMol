import React from "react";
import Box from "@mui/material/Box";
import Button from "@mui/material/Button";
import Card from "@mui/material/Card";
import CardContent from "@mui/material/CardContent";
import Stack from "@mui/material/Stack";
import Typography from "@mui/material/Typography";
import MuiLink from "@mui/material/Link";
import Tooltip from "@mui/material/Tooltip";
import NotificationsRoundedIcon from "@mui/icons-material/NotificationsRounded";
import RefreshIcon from "@mui/icons-material/Refresh";
import { useTheme } from "@mui/material/styles";
import { useNotifications } from "../NotificationProvider";
import { Link as RouterLink } from "react-router-dom";
import { DialogImportCompound } from "./DialogImportCompound";
import { DialogImportGeneCluster } from "./DialogImportGeneClusters";
import { WorkspaceItemCard } from "./WorkspaceItemCard";
import { Session } from "../../features/session/types";
import { deleteSessionItem, refreshSession } from "../../features/session/api";
import { NewCompoundJob } from "../../features/jobs/types";
import { MAX_ITEMS, importCompound, importCompoundsBatch, importClustersBatch } from "../../features/jobs/api";

const MAX_FILE_SIZE_MB = 2;
const MAX_FILE_SIZE_BYTES = MAX_FILE_SIZE_MB * 1024 * 1024;

async function parseCompoundFile(file: File, matchStereochemistry: boolean): Promise<NewCompoundJob[]> {
  const text = await file.text();

  // Normalize newlines and split lines
  const lines = text.replace(/\r\n/g, "\n").replace(/\r/g, "\n").split("\n");
  if (lines.length === 0) return [];

  const headerLine = lines[0].trim();

  // Detect delimiter based on file extension
  let delimiter = ",";
  if (file.name.endsWith(".tsv") || file.name.endsWith(".txt")) {
    delimiter = "\t";
  } else if (file.name.endsWith(".csv")) {
    delimiter = ",";
  };

  const headers = headerLine.split(delimiter).map(h => h.trim().toLowerCase());
  const nameIdx = headers.indexOf("name");
  const smilesIdx = headers.indexOf("smiles");

  if (nameIdx === -1 || smilesIdx === -1) {
    throw new Error("File must contain 'name' and 'smiles' columns in the header.");
  };

  const compounds: NewCompoundJob[] = [];

  for (let i = 1; i < lines.length; i++) {
    const line = lines[i].trim();
    if (line === "") continue; // skip empty lines

    const cols = line.split(delimiter).map(c => c.trim());
    const name = (cols[nameIdx] ?? "").trim();
    const smiles = (cols[smilesIdx] ?? "").trim();

    if (!name || !smiles) continue;

    compounds.push({ name, smiles, matchStereochemistry });
  };

  return compounds;
};

type WorkspaceUploadProps = {
  session: Session;
  setSession: React.Dispatch<React.SetStateAction<Session | null>>;
};

export const WorkspaceUpload: React.FC<WorkspaceUploadProps> = ({ session, setSession }) => {
  const theme = useTheme();
  const { pushNotification } = useNotifications();

  const [openCompounds, setOpenCompounds] = React.useState(false);
  const [openClusters, setOpenClusters] = React.useState(false);

  const [selectedIds, setSelectedIds] = React.useState<Set<string>>(new Set());
  const [deletingIds, setDeletingIds] = React.useState<Set<string>>(new Set());
  const [refreshSpinKey, setRefreshSpinKey] = React.useState(0); // to force re-mount of refresh icon

  // Clean up deletingIds when session items change
  React.useEffect(() => {
    setDeletingIds((prev) => {
      const liveIds = new Set(session.items.map(it => it.id));
      const next = new Set<string>();

      prev.forEach((id) => {
        if (liveIds.has(id)) next.add(id);
      });

      return next;
    });
  }, [session.items]);

  // Wrap parent setter (Session | null) into the deps shape (Session-only functional updater)
  const setSessionSafe = React.useCallback(
    (updater: (prev: Session) => Session) => {
      setSession((prev) => (prev ? updater(prev) : prev));
    },
    [setSession]
  );

  // Helper to build deps for import service
  const deps = React.useMemo(
    () => ({
      setSession: setSessionSafe,
      pushNotification,
      sessionId: session.sessionId,
    }),
    [setSessionSafe, pushNotification, session.sessionId]
  );

  // Helper for manual refresh
  const handleManualRefresh = async () => {
    setRefreshSpinKey((k) => k + 1); // trigger re-mount of icon to restart animation
    try {
      const fresh = await refreshSession(session.sessionId);
      setSession(() => fresh);
      pushNotification("Workspace session refreshed successfully!", "success");
    } catch (err) {
      const msg = err instanceof Error ? err.message : String(err);
      pushNotification(`Failed to refresh session: ${msg}`, "error");
    };
  };

  // Selection helpers
  const toggleSelectItem = (id: string) => {
    // Prevent toggling if deleting
    if (deletingIds.has(id)) return;

    setSelectedIds(prev => {
      const next = new Set(prev);
      if (next.has(id)) next.delete(id);
      else next.add(id);
      return next;
    })
  };

  const handleSelectAll = () => {
    if (!session.items.length) return;
    setSelectedIds(new Set(session.items.map(item => item.id)));
  };

  const handleClearSelection = () => {
    setSelectedIds(new Set());
  };

  const handleDeleteItem = async (id: string) => {
    // Mark as deleting (UI only)
    setDeletingIds(prev => new Set(prev).add(id));

    try {
      await deleteSessionItem(session.sessionId, id);
      // DO NOTHING ELSE
      // SSE refresh will remove item from session
    } catch (err) {
      const msg = err instanceof Error ? err.message : String(err);
      pushNotification(`Failed to delete item: ${msg}`, "error");
      setDeletingIds(prev => {
        const n = new Set(prev);
        n.delete(id);
        return n;
      });
    }
  };

  const handleDeleteSelected = async () => {
    if (selectedIds.size === 0) return;

    const ids = Array.from(selectedIds);
    setSelectedIds(new Set()); // UI-only

    try {
      for (const id of ids) {
        await deleteSessionItem(session.sessionId, id);
      }
      // SSE will update the session state
    } catch (err) {
      const msg = err instanceof Error ? err.message : String(err);
      pushNotification(`Failed to delete selected items from session: ${msg}`, "error");
    };
  };

  // Open dailogs
  const handleOpenCompounds = (event: React.MouseEvent<HTMLButtonElement, MouseEvent>) => {
    event.currentTarget.blur(); // prevents 'Blocked aria-hidden on an element' warning
    setOpenCompounds(true);
  };

  const handleOpenBGCs = (event: React.MouseEvent<HTMLButtonElement, MouseEvent>) => {
    event.currentTarget.blur(); // prevents 'Blocked aria-hidden on an element' warning
    setOpenClusters(true);
  };

  // Import compound handlers
  const handleImportSingleCompound = async({ name, smiles, matchStereochemistry}: { name: string; smiles: string; matchStereochemistry: boolean }) => {
    await importCompound(deps, { name, smiles, matchStereochemistry });
  };

  const handleImportBatchCompounds = async (file: File, matchStereochemistry: boolean) => {
    if (file.size > MAX_FILE_SIZE_BYTES) {
      pushNotification(`The file "${file.name}" exceeds the maximum size of ${MAX_FILE_SIZE_MB} MB and was not imported.`, "error");
      return;
    };

    try {
      const compounds = await parseCompoundFile(file, matchStereochemistry);
      await importCompoundsBatch(deps, compounds);
    } catch (err) {
      const msg = err instanceof Error ? err.message : String(err);
      pushNotification(`Failed to parse compound file: ${msg}`, "error");
    };
  };

  // Import cluster handler
  const handleImportClusters = async (files: File[]) => {
    if (!files.length) return;

    const oversized = files.filter(f => f.size > MAX_FILE_SIZE_BYTES);
    if (oversized.length > 0) {
      pushNotification(`Some files exceed the maximum size of ${MAX_FILE_SIZE_MB} MB and were not imported: ${oversized.map(f => f.name).join(", ")}`, "error");

      // Keep only files within size limit
      files = files.filter(f => f.size <= MAX_FILE_SIZE_BYTES);
    };

    // Check if any files remain after filtering on file size
    if (files.length === 0) {
      pushNotification("No valid files to import after size filtering.", "warning");
      return;
    };

    let payloads: { name: string; fileContent: string }[] = [];

    try {
      payloads = await Promise.all(
        files.map(async (file) => ({
          name: file.name,
          fileContent: await file.text(),
        }))
      )
      await importClustersBatch(deps, payloads);
    } catch (err) {
      const msg = err instanceof Error ? err.message : String(err);
      pushNotification(`Failed to import BGC files: ${msg}`, "error");
    };
  };

  // Selection states
  const anySelected = selectedIds.size > 0;
  const allSelected = session.items.length > 0 && selectedIds.size === session.items.length;

  return (
    <Box
      sx={{
          width: "100%",
          mx: "auto",
          display: "flex",
          flexDirection: "column",
          gap: "16px"
      }}
    >
      <Card
        variant="outlined"
        sx={{
          display: "flex",
          flexDirection: "column",
          gap: "8px",
          flexGrow: 1,
        }}
      >
        <CardContent>
          <Typography component="h1" variant="subtitle1">
            Getting started
          </Typography>
          <Typography variant="body1">
            In this tab you can import compounds and biosynthetic gene clusters (BGCs) into your workspace. Use the buttons below to upload your data files. After importing, you can visualize and analyze your data within the&nbsp;
            <MuiLink
              component={RouterLink}
              to="/dashboard/discovery"
              underline="hover"
              color={(theme.vars || theme).palette.primary.main}
              sx={{ fontWeight: "500" }}
            >
              Discovery
            </MuiLink>
            &nbsp;and&nbsp;
            <MuiLink
              component={RouterLink}
              to="/dashboard/enrichment"
              underline="hover"
              color={(theme.vars || theme).palette.primary.main}
              sx={{ fontWeight: "500" }}
            >
              Enrichment
            </MuiLink>
            &nbsp;tabs. A maximum of <b>{MAX_ITEMS} items</b> can be imported into the workspace. Keep an eye on <NotificationsRoundedIcon fontSize={'small'} sx={{ verticalAlign: 'middle' }} /> for updates on your queries.
          </Typography>

          <Stack direction="row" spacing={2} sx={{ mt: 2 }}>
            <Button variant="contained" onClick={handleOpenCompounds}>
              Import compounds
            </Button>

            <Button variant="contained" onClick={handleOpenBGCs} disabled>
              Import BGCs
            </Button>
          </Stack>
        </CardContent>
      </Card>

      <DialogImportCompound
        open={openCompounds}
        onClose={() => setOpenCompounds(false)}
        onImportSingle={handleImportSingleCompound}
        onImportBatch={handleImportBatchCompounds}
      />

      <DialogImportGeneCluster
        open={openClusters}
        onClose={() => setOpenClusters(false)}
        onImport={handleImportClusters}
      />

      {session.items.length > 0 && (
        <Card variant="outlined">
          <CardContent>
            <Stack
              direction="row"
              alignItems="center"
              justifyContent="space-between"
              sx={{ mb: 1.5 }}
            >
              <Stack direction="row" alignItems="center" spacing={1}>
                <Typography component="h1" variant="subtitle1">
                  Workspace items ({session.items.length}/{MAX_ITEMS})
                </Typography>
                <Tooltip title="Refresh workspace items" arrow>
                  <RefreshIcon
                    key={refreshSpinKey} // to force re-mount
                    fontSize="small"
                    sx={{
                      "@keyframes refresh-spin": {
                        from: { transform: "rotate(0deg)" },
                        to: { transform: "rotate(360deg)" },
                      },
                      animation: refreshSpinKey > 0 ? "refresh-spin 0.6s linear" : "none",
                      cursor: "pointer",
                      color: (theme.vars || theme).palette.text.secondary,
                      "&:hover": {
                        color: (theme.vars || theme).palette.text.primary,
                      },
                    }}
                    onClick={handleManualRefresh}
                  />
                </Tooltip>
              </Stack>
              <Stack direction="row" spacing={1}>
                <Button
                  size="small"
                  variant="text"
                  onClick={handleSelectAll}
                  disabled={session.items.length === 0 || allSelected}
                >
                  Select all
                </Button>
                <Button
                  size="small"
                  variant="text"
                  onClick={handleClearSelection}
                  disabled={!anySelected}
                >
                  Clear selection
                </Button>
                <Button
                  size="small"
                  variant="text"
                  color="error"
                  onClick={handleDeleteSelected}
                  disabled={!anySelected}
                >
                  Delete selected
                </Button>
              </Stack>
            </Stack>

            <Stack spacing={1}>
              {session.items.map((item) => (
                <WorkspaceItemCard
                  key={item.id}
                  sessionId={session.sessionId}
                  item={item}
                  selected={selectedIds.has(item.id)}
                  disabled={deletingIds.has(item.id)}
                  onToggleSelect={toggleSelectItem}
                  onDelete={handleDeleteItem}
                />
              ))}
            </Stack>
          </CardContent>
        </Card>
      )}

    </Box>
  );
};