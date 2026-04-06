import React from "react";
import Box from '@mui/material/Box';
import Button from '@mui/material/Button';
import Card from '@mui/material/Card';
import CardContent from '@mui/material/CardContent';
import Stack from '@mui/material/Stack';
import Typography from '@mui/material/Typography';
import MuiLink from "@mui/material/Link";
import NotificationsRoundedIcon from '@mui/icons-material/NotificationsRounded';
import { useTheme } from "@mui/material/styles";
import { useNotifications } from "../NotificationProvider";
import { Link as RouterLink } from "react-router-dom";
import { DialogImportCompound } from "./DialogImportCompound";
import { DialogImportGeneCluster } from "./DialogImportGeneClusters"
import { DialogViewItem } from "./DialogViewItem";
import { WorkspaceItemCard } from "./WorkspaceItemCard";
import { Session } from "../../features/session/types";
import { NewCompoundJob } from "../../features/jobs/types";
import {
  MAX_ITEMS,
  importCompound,
  importCompoundsBatch,
  importGeneClustersBatch,
} from "../../features/jobs/api";

// const MAX_ITEMS = 200;
const MAX_FILE_SIZE_MB = 2;
const MAX_FILE_SIZE_BYTES = MAX_FILE_SIZE_MB * 1024 * 1024;

async function parseCompoundFile(file: File): Promise<NewCompoundJob[]> {
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
  }

  const headers = headerLine.split(delimiter).map(h => h.trim().toLowerCase());
  const nameIdx = headers.indexOf("name");
  const smilesIdx = headers.indexOf("smiles");

  if (nameIdx === -1 || smilesIdx === -1) {
    throw new Error("File must contain 'name' and 'smiles' columns in the header.");
  }

  const compounds: NewCompoundJob[] = [];

  for (let i = 1; i < lines.length; i++) {
    const line = lines[i].trim();
    if (line === "") continue; // skip empty lines

    const cols = line.split(delimiter).map(c => c.trim());
    const name = (cols[nameIdx] ?? "").trim();
    const smiles = (cols[smilesIdx] ?? "").trim();

    if (!name || !smiles) continue;

    compounds.push({ name, smiles });
  }

  return compounds;
}

type WorkspaceUploadProps = {
  session: Session;
  setSession: (updated: (prev: Session) => Session) => void;
}

export const WorkspaceUpload: React.FC<WorkspaceUploadProps> = ({ session, setSession }) => {
  const theme = useTheme();
  const { pushNotification } = useNotifications();

  const [openCompounds, setOpenCompounds] = React.useState(false);
  const [openGeneClusters, setOpenGeneClusters] = React.useState(false);
  const [openView, setOpenView] = React.useState(false);
  const [viewingItemId, setViewingItemId] = React.useState<string | null>(null);

  const [selectedIds, setSelectedIds] = React.useState<Set<string>>(new Set());

  // Helper to build deps for import service
  const deps = React.useMemo(
    () => ({
      setSession,
      pushNotification,
      sessionId: session.sessionId,
    }),
    [setSession, session.sessionId]
  )

  // Renaming helper
  const handleRenameItem = (id: string, newName: string) => {
    setSession((prev) => ({
      ...prev,
      items: prev.items.map((item) =>
        item.id === id
          ? { ...item, name: newName }
          : item
      ),
    }))
  }

  // Viewing helper
  const handleViewItem = (id: string) => {
    setViewingItemId(id);
    setOpenView(true);
  }

  // Selection helpers
  const toggleSelectItem = (id: string) => {
    setSelectedIds(prev => {
      const next = new Set(prev);
      if (next.has(id)) next.delete(id);
      else next.add(id);
      return next;
    })
  }

  const handleDeleteItem = (id: string) => {
    setSession(prev => ({
      ...prev,
      items: prev.items.filter(item => item.id !== id),
    }))
    setSelectedIds(prev => {
      const next = new Set(prev);
      next.delete(id);
      return next;
    })
  }

  const handleSelectAll = () => {
    if (!session.items.length) return;
    setSelectedIds(new Set(session.items.map(item => item.id)));
  }

  const handleClearSelection = () => {
    setSelectedIds(new Set());
  }

  const handleDeleteSelected = () => {
    if (selectedIds.size === 0) return;
    const toDelete = new Set(selectedIds);

    setSession(prev => ({
      ...prev,
      items: prev.items.filter(item => !toDelete.has(item.id)),
    }))
    setSelectedIds(new Set());
  }

  const handleImportSingleCompound = async({ name, smiles}: { name: string; smiles: string }) => {
    await importCompound(deps, { name, smiles });
  }

  const handleImportBatchCompounds = async (file: File) => {
    if (file.size > MAX_FILE_SIZE_BYTES) {
      pushNotification(`The file "${file.name}" exceeds the maximum size of ${MAX_FILE_SIZE_MB} MB and was not imported.`, "error");
      return;
    }

    try {
      const compounds = await parseCompoundFile(file);
      await importCompoundsBatch(deps, compounds);
    } catch (err) {
      const msg = err instanceof Error ? err.message : String(err);
      pushNotification(`Failed to parse compound file: ${msg}`, "error");
    }
  }

  const handleImportGeneClusters = async (files: File[]) => {
    if (!files.length) return;

    const oversized = files.filter(f => f.size > MAX_FILE_SIZE_BYTES);
    if (oversized.length > 0) {
      pushNotification(`The following files exceed the maximum size of ${MAX_FILE_SIZE_MB} MB and were not imported: ${oversized.map(f => f.name).join(", ")}`, "error");

      // Keep only files within size limit
      files = files.filter(f => f.size <= MAX_FILE_SIZE_BYTES);
    }

    // Check if any files remain after filtering on file size
    if (files.length === 0) {
      pushNotification("No gene cluster files to import after size check.", "warning");
      return;
    }

    let payloads: { name: string; fileContent: string }[] = [];
    try {
      payloads = await Promise.all(
        files.map(async (file) => ({
          name: file.name,
          fileContent: await file.text(),
        }))
      )
    } catch (err) {
      const msg = err instanceof Error ? err.message : String(err);
      pushNotification(`Failed to read gene cluster files: ${msg}`, "error");
      return;
    }

    await importGeneClustersBatch(deps, payloads);
  }

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
            To get started, you can import compounds and gene clusters into your workspace. Use the buttons below to upload your data files. After importing, you can visualize and analyze your data within the&nbsp;
            <MuiLink
              component={RouterLink}
              to="/dashboard/explore"
              underline="hover"
              color={(theme.vars || theme).palette.primary.main}
              sx={{ fontWeight: "500" }}
            >
              Explore tab
            </MuiLink>
            . A maximum of <b>{MAX_ITEMS} items</b> can be imported into the workspace. Keep an eye on <NotificationsRoundedIcon fontSize={'small'} sx={{ verticalAlign: 'middle' }} /> for updates on your queries.
          </Typography>

          <Stack direction="row" spacing={2} sx={{ mt: 2 }}>
            <Button variant="contained" onClick={() => setOpenCompounds(true)}>
              Import compounds
            </Button>
            <Button variant="contained" onClick={() => setOpenGeneClusters(true)}>
              Import gene clusters
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
        open={openGeneClusters}
        onClose={() => setOpenGeneClusters(false)}
        onImport={handleImportGeneClusters}
      />

      <DialogViewItem
        open={openView}
        item={session.items.find(item => item.id === viewingItemId) ?? null}
        onClose={() => {
          setOpenView(false);
          setViewingItemId(null);
        }}
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
              <Typography component="h1" variant="subtitle1">
                Workspace items ({session.items.length}/{MAX_ITEMS})
              </Typography>
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
                  item={item}
                  selected={selectedIds.has(item.id)}
                  onToggleSelect={toggleSelectItem}
                  onDelete={handleDeleteItem}
                  onView={handleViewItem}
                  onRename={handleRenameItem}
                />
              ))}
            </Stack>
          </CardContent>
        </Card>
      )}

    </Box>
  )
}