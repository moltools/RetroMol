import React from "react";
import Alert from "@mui/material/Alert";
import CircularProgress from "@mui/material/CircularProgress";
import Box from "@mui/material/Box";
import Stack from "@mui/material/Stack";
import Tooltip from "@mui/material/Tooltip";
import Typography from "@mui/material/Typography";
import InfoOutlinedIcon from "@mui/icons-material/InfoOutlined";
import { useQuery } from "@tanstack/react-query";
import { Session, SessionItem } from "../../features/session/types";
import { reconstructCompound } from "../../features/reconstruction/api";
import { getClusterReadout } from "../../features/clusters/api";
import { DialogWindow } from "../DialogWindow";
import { ErrorBoundary } from "../ErrorBoundary";
import SmilesDrawerContainer from "../SmilesDrawerContainer.js";
import { PrimarySequenceRows, usePrimarySequenceEditor } from "./PrimarySequenceEditor";
import { ClusterReadoutRows } from "./ClusterReadoutRows";

type HighlightAtom = [number, string];

type DialogViewItemProps = {
  session: Session;
  setSession: React.Dispatch<React.SetStateAction<Session | null>>;
  item: SessionItem;
  open: boolean;
  onClose: () => void;
};

function AnnotatedArrow({ annotation }: { annotation: string }) {
  return (
    <Box
      sx={{
        display: 'flex',
        flexDirection: 'column',
        alignItems: 'center',
        justifyContent: 'center',
        flexShrink: 0, // don’t let it compress
      }}
    >
      <Typography variant='caption' gutterBottom>
        {annotation}
      </Typography>
      <Box
        component='svg'
        width={140} // total width of arrow+shaft
        height={32}
        viewBox='0 0 140 32'
      >
        {/* horizontal line */}
        <line
          x1='0'
          y1='16'
          x2='120'
          y2='16'
          stroke='currentColor'
          strokeWidth='2'
        />
        {/* arrow head */}
        <polygon
          points='120,8 120,24 140,16'
          fill='currentColor'
        />
      </Box>
    </Box>
  );
};

function DescriptionBox({ title, description }: { title: string; description: string }) {
  return (
    <Box sx={{ width: '100%' }}>
      <Typography variant='subtitle2' fontWeight='bold'>
        {title}
      </Typography>
      <Typography variant='body2'>
        {description}
      </Typography>
    </Box>
  );
};

export const DialogViewItem: React.FC<DialogViewItemProps> = ({
  session,
  setSession,
  item,
  open,
  onClose,
}) => {
  const sessionId = session.sessionId;
  const isCompound = item.kind === "compound"; // there are only two types: "compound" and "cluster"
  const itemScore = typeof item.score === "number" ? item.score : 0;

  const [selectedTags, setSelectedTags] = React.useState<number[]>([]);

  const handleToggleMotif = (tags: number[]) => {
    setSelectedTags((prev) => {
      const allSelected = tags.every((tag) => prev.includes(tag));
      if (allSelected) {
        return prev.filter((tag) => !tags.includes(tag));
      }
      return Array.from(new Set([...prev, ...tags]));
    })
  }

  // Memoized so it's a stable reference across re-renders that don't touch
  // selectedTags — otherwise SmilesDrawerContainer redraws on every unrelated
  // parent re-render (e.g. the item card's periodic "updated Xs ago" tick).
  const highlightAtoms = React.useMemo<HighlightAtom[]>(
    () => selectedTags.map((tag) => [tag, "#027bf3"]),
    [selectedTags]
  );

  // Clear selection when dialog/item changes
  React.useEffect(() => {
    if (open) {
      setSelectedTags([]);
    }
  }, [open, item.id]);

  const reconstructionQuery = useQuery({
    queryKey: ["reconstructCompound", sessionId, item.id],
    queryFn: ({ signal }) => reconstructCompound(sessionId, item.id, signal),
    enabled: open && isCompound,
  });

  const data = reconstructionQuery.data ?? null;
  const loading = reconstructionQuery.isLoading;
  const error = reconstructionQuery.error
    ? (reconstructionQuery.error as Error).message || "Unknown error"
    : null;

  const clusterReadoutQuery = useQuery({
    queryKey: ["getClusterReadout", sessionId, item.id],
    queryFn: ({ signal }) => getClusterReadout(sessionId, item.id, signal),
    enabled: open && !isCompound && item.status === "done",
  });

  // resetSignal is `open` -- re-seeding on every open (even for the same item)
  // matches the dialog's existing "fresh state each time" behavior for selectedTags.
  const editor = usePrimarySequenceEditor(session, setSession, item, data, open);
  const { editing, setEditing, anyDirty, saving, handleCancelEdit, handleSaveAll } = editor;

  return (
    <DialogWindow
      open={open}
      onClose={onClose}
      title="View item"
      dividers
      dirty={editing && anyDirty}
      actions={[
        ...(isCompound && editing
          ? [
              { key: "cancel-edit", label: "Cancel", variant: "text" as const, color: "inherit" as const, onClick: handleCancelEdit },
              {
                key: "save",
                label: saving ? "Saving..." : "Save changes",
                variant: "contained" as const,
                color: "primary" as const,
                onClick: handleSaveAll,
                disabled: !anyDirty || saving,
                startIcon: saving ? <CircularProgress size={16} color="inherit" /> : undefined,
              },
            ]
          : []),
        ...(isCompound && !editing
          ? [
              {
                key: "edit",
                label: "Edit sequences",
                variant: "outlined" as const,
                color: "primary" as const,
                onClick: () => setEditing(true),
                disabled: loading || !data || data.length === 0,
              },
            ]
          : []),
        { key: "close", label: "Close", variant: "text" as const, color: "inherit" as const, onClick: onClose },
      ]}
      maxWidth={"lg"}
    >
      {!isCompound && item.kind === "cluster" && (
        <>
          {item.status === "queued" && (
            <Typography variant="body2" color="text.secondary">
              Waiting to be parsed...
            </Typography>
          )}

          {item.status === "processing" && <CircularProgress size={24} />}

          {item.status === "error" && (
            <Alert severity="error">
              {item.errorMessage || "Parsing failed."}
            </Alert>
          )}

          {item.status === "done" && clusterReadoutQuery.isLoading && <CircularProgress size={24} />}

          {item.status === "done" && clusterReadoutQuery.error && (
            <Alert severity="error">
              {(clusterReadoutQuery.error as Error).message || "Failed to load parsed gene cluster."}
            </Alert>
          )}

          {item.status === "done" && clusterReadoutQuery.data && (
            <ClusterReadoutRows readouts={clusterReadoutQuery.data.readouts} />
          )}
        </>
      )}

      {loading && (
        <CircularProgress size={24} />
      )}

      {error && (
        <Alert severity="error" sx={{ mb: 2 }}>
          {error}
        </Alert>
      )}

      {(isCompound && !loading && !error) && (
        <Box sx={{ display: "flex", flexDirection: "column", gap: 2 }}>
          {(itemScore < 0.5 || (data?.length ?? 0) === 0) && (
            <Stack direction="row" spacing={0.5} alignItems="center">
              <Tooltip
                title="Coverage is the share of the molecule RetroMol could match to known building blocks. Low or 0% coverage means most (or all) of the structure falls outside its rule set, e.g. an unusual scaffold or a modification the current rules don't recognize. A low coverage leads to primary sequences that may be sparse, incomplete, or empty."
                placement="bottom-start"
                arrow
              >
                <Stack direction="row" spacing={0.5} alignItems="center" sx={{ cursor: "help", width: "fit-content" }}>
                  <InfoOutlinedIcon fontSize="small" color="warning" />
                  <Typography variant="caption" color="text.secondary">
                    Low coverage -- parts of the structure below may be sparse or empty
                  </Typography>
                </Stack>
              </Tooltip>
            </Stack>
          )}
          <Box
            sx={{
              display: 'flex',
              flexDirection: 'row',
              gap: 2,
              alignItems: 'center',
              overflowX: 'auto',
              overflowY: 'hidden',
              scrollbarGutter: 'stable',
              '&::-webkit-scrollbar': {
                height: 12,
              },
              '&::-webkit-scrollbar-track': {
                background: '#f1f1f1',
                borderRadius: 4,
              },
              '&::-webkit-scrollbar-thumb': {
                backgroundColor: '#ccc',
                borderRadius: 4,
                '&:hover': {
                  backgroundColor: '#aaa',
                  cursor: 'pointer',
                },
              },
            }}
          >
            <Box
              sx={{
                flex: 1,
                display: 'flex',
                justifyContent: 'center',
                alignItems: 'center',
                pl: 2
              }}
            >
              <ErrorBoundary what="molecule structure" fallback={
                <Typography variant="caption" color="text.secondary">Could not render this structure.</Typography>
              }>
                <SmilesDrawerContainer
                  identifier={`smiles-drawer-${sessionId}-${item.id}-full`}
                  smiles={data?.[0]?.tagged_input_smiles ?? ""}
                  size={300}
                  highlightAtoms={highlightAtoms}
                />
              </ErrorBoundary>
            </Box>
            <AnnotatedArrow annotation={'Linearization'} />
            <Box
              sx={{
                flex: 1,
                display: 'flex',
                justifyContent: 'center',
                alignItems: 'center'
              }}
            >
              <Box
                sx={{
                  flex: 1,
                  display: "flex",
                  flexDirection: "column",
                  gap: 1,
                  justifyContent: "center",
                  alignItems: "center",
                }}
              >
                {(data ?? []).map((reconstruction, idx) => (
                  <React.Fragment key={`backbone-fragment-${idx}`}>
                    <ErrorBoundary what="molecule structure" fallback={
                      <Typography variant="caption" color="text.secondary">Could not render this structure.</Typography>
                    }>
                      <SmilesDrawerContainer
                        identifier={`smiles-drawer-${sessionId}-${item.id}-preprocessed-${idx}`}
                        smiles={reconstruction.tagged_backbone_smiles}
                        size={260}
                        highlightAtoms={highlightAtoms}
                      />
                    </ErrorBoundary>

                    {idx < (data?.length ?? 0) - 1 && (
                      <Typography
                        variant="h6"
                        sx={{
                          lineHeight: 1,
                          color: "text.secondary",
                          fontWeight: 500,
                        }}
                      >
                        +
                      </Typography>
                    )}
                  </React.Fragment>
                ))}
              </Box>
            </Box>
            <AnnotatedArrow annotation={'Sequencing'} />
            <Box
              sx={{
                // Sized to its natural content width (not squeezed by the other
                // flex:1 siblings) -- the full sequence should be readable without
                // an extra inner scrollbar; the outer row above already scrolls
                // horizontally for anything that doesn't fit.
                flex: '0 0 auto',
                display: 'flex',
                justifyContent: 'center',
                alignItems: 'center',
                pr: 2
              }}
            >
              <Box
                sx={{
                  flex: 1,
                  display: "flex",
                  flexDirection: "column",
                  gap: 2,
                  justifyContent: "center",
                  alignItems: "center",
                  pr: 2,
                }}
              >
                <Box
                  sx={{
                    flex: 1,
                    display: "flex",
                    flexDirection: "column",
                    gap: 2,
                    alignItems: "flex-start",
                    pr: 2,
                  }}
                >
                  <PrimarySequenceRows
                    item={item}
                    data={data ?? []}
                    state={editor}
                    selectedTags={selectedTags}
                    onToggleMotif={handleToggleMotif}
                  />
                </Box>
              </Box>
            </Box>
          </Box>
          <DescriptionBox
            title={'Explanation'}
            description={'The input SMILES structure (left) is processed by RetroMol into non-overlapping building blocks, and from these building blocks a linear backbone is reconstructed (middle). The primary sequence, a text representation of the linear backbone, is seen on the right. You can highlight individual motifs by clicking them in the primary sequence above. If a sequence was parsed wrong or incompletely, use "Edit sequences" to fix it by hand. Saved edits are kept with this compound and reused when you query it from the Discovery tab; "Revert" on a row discards the edit and restores the algorithm’s own parse.'}
          />
        </Box>
      )}
    </DialogWindow>
  );
};