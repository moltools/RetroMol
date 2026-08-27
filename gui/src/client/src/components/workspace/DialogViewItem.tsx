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
import { reconstructGeneCluster } from "../../features/clusters/api";
import { DialogWindow } from "../DialogWindow";
import { ErrorBoundary } from "../ErrorBoundary";
import { ExportImageButton } from "../ExportImageButton";
import SmilesDrawerContainer from "../SmilesDrawerContainer.js";
import { DrawingAttribution } from "../DrawingAttribution";
import { PrimarySequenceOverview, PrimarySequenceRows, usePrimarySequenceEditor } from "./PrimarySequenceEditor";
import { AnnotatedArrow } from "./AnnotatedArrow";
import { ClusterReadoutDiagram, useClusterPrimarySequenceEditor } from "./ClusterReadoutDiagram";

type HighlightAtom = [number, string];

type DialogViewItemProps = {
  session: Session;
  setSession: React.Dispatch<React.SetStateAction<Session | null>>;
  item: SessionItem;
  open: boolean;
  onClose: () => void;
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
  const diagramRef = React.useRef<HTMLDivElement>(null);

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

  const data = reconstructionQuery.data?.reconstructions ?? null;
  const dbMatchingSequences = reconstructionQuery.data?.dbMatchingSequences ?? [];
  const loading = reconstructionQuery.isLoading;
  const error = reconstructionQuery.error
    ? (reconstructionQuery.error as Error).message || "Unknown error"
    : null;

  // Four possible outcomes from the backend (see retromol_synthesis.reconstruction):
  //  1. hasReconstructions && allHaveBackbone: full structure -> backbone -> sequence flow.
  //  2. hasReconstructions && !isUnordered && !allHaveBackbone: a primary sequence was
  //     found, but the hardcoded backbone-fusion chemistry couldn't rebuild a backbone
  //     structure for it -- skip the backbone step and go straight from structure to
  //     sequence.
  //  3. hasReconstructions && isUnordered: building blocks were identified individually
  //     (e.g. a branched or cyclic assembly), but couldn't be threaded into a single
  //     order -- show them as an unordered set, not a sequence.
  //  4. !hasReconstructions: nothing identified at all -- show the structure on its own.
  const hasReconstructions = (data?.length ?? 0) > 0;
  const isUnordered = hasReconstructions && (data ?? []).every((r) => r.ordered === false);
  const allHaveBackbone = hasReconstructions && !isUnordered && (data ?? []).every((r) => !!r.tagged_backbone_smiles);
  const backboneWarning = (data ?? []).find((r) => r.backbone_warning)?.backbone_warning ?? null;

  const clusterReconstructionQuery = useQuery({
    queryKey: ["reconstructGeneCluster", sessionId, item.id],
    queryFn: ({ signal }) => reconstructGeneCluster(sessionId, item.id, signal),
    enabled: open && !isCompound && item.status === "done",
  });
  const clusterData = clusterReconstructionQuery.data ?? null;

  // resetSignal is `open` -- re-seeding on every open (even for the same item)
  // matches the dialog's existing "fresh state each time" behavior for selectedTags.
  const editor = usePrimarySequenceEditor(session, setSession, item, data, open);
  const { editing, setEditing, anyDirty, saving, handleCancelEdit, handleSaveAll } = editor;

  const clusterEditor = useClusterPrimarySequenceEditor(session, setSession, item, clusterData, open);
  const {
    editing: clusterEditing,
    setEditing: setClusterEditing,
    anyDirty: clusterAnyDirty,
    saving: clusterSaving,
    handleCancelEdit: handleClusterCancelEdit,
    handleSaveAll: handleClusterSaveAll,
  } = clusterEditor;
  const clusterDiagramRef = React.useRef<HTMLDivElement>(null);

  return (
    <DialogWindow
      open={open}
      onClose={onClose}
      title="View item"
      dividers
      dirty={(editing && anyDirty) || (clusterEditing && clusterAnyDirty)}
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
                variant: "contained" as const,
                color: "primary" as const,
                onClick: () => setEditing(true),
                disabled: loading || !data || data.length === 0 || isUnordered,
              },
            ]
          : []),
        ...(item.kind === "cluster" && clusterEditing
          ? [
              { key: "cancel-cluster-edit", label: "Cancel", variant: "text" as const, color: "inherit" as const, onClick: handleClusterCancelEdit },
              {
                key: "save-cluster",
                label: clusterSaving ? "Saving..." : "Save changes",
                variant: "contained" as const,
                color: "primary" as const,
                onClick: handleClusterSaveAll,
                disabled: !clusterAnyDirty || clusterSaving,
                startIcon: clusterSaving ? <CircularProgress size={16} color="inherit" /> : undefined,
              },
            ]
          : []),
        ...(item.kind === "cluster" && !clusterEditing
          ? [
              {
                key: "edit-cluster",
                label: "Edit sequences",
                variant: "contained" as const,
                color: "primary" as const,
                onClick: () => setClusterEditing(true),
                disabled: clusterReconstructionQuery.isLoading || !clusterData || clusterData.length === 0,
              },
            ]
          : []),
        { key: "close", label: "Close", variant: "text" as const, color: "inherit" as const, onClick: onClose },
      ]}
      maxWidth={"lg"}
    >
      {!isCompound && item.kind === "cluster" && (
        <Box sx={{ display: "flex", flexDirection: "column", gap: 2 }}>
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

          {item.status === "done" && clusterReconstructionQuery.isLoading && <CircularProgress size={24} />}

          {item.status === "done" && clusterReconstructionQuery.error && (
            <Alert severity="error">
              {(clusterReconstructionQuery.error as Error).message || "Failed to load parsed gene cluster."}
            </Alert>
          )}

          {item.status === "done" && clusterData && (
            <>
              <Box sx={{ display: "flex", justifyContent: "flex-end" }}>
                <ExportImageButton
                  targetRef={clusterDiagramRef}
                  filename={`retromol-${(item.name || item.id).replace(/[^a-z0-9]+/gi, "-")}-readout`}
                  label="Download the diagram below as a PNG"
                />
              </Box>
              <Box ref={clusterDiagramRef} sx={{ p: 1 }}>
                <ClusterReadoutDiagram item={item} data={clusterData} state={clusterEditor} />
              </Box>
              <DescriptionBox
                title="Explanation"
                description='Each antiSMASH region is shown as a "genes & domains" row (left), the ordered NRPS/PKS modules RetroMol detected, and a "readout" row (right) -- the same primary-sequence vocabulary a compound is parsed into. Click a chip on either row to highlight the module(s)/token it corresponds to on the other. If a readout was resolved wrong or incompletely, use "Edit sequences" to fix it by hand; saved edits are kept with this gene cluster and reused when you query it from the Discovery tab. "Revert" on a row discards the edit and restores the algorithm’s own parse.'
              />
            </>
          )}
        </Box>
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
          {!hasReconstructions && (
            <Alert severity="warning">
              This compound couldn't be parsed with RetroMol's current rule set, so no
              primary sequence could be derived from it. The structure is shown below as-is.
            </Alert>
          )}

          {hasReconstructions && !allHaveBackbone && (
            <Alert severity="warning">
              {backboneWarning ??
                (isUnordered
                  ? "RetroMol identified these building blocks individually, but couldn't determine a biosynthetic order between them."
                  : "The linear backbone could not be reconstructed for this structure.")}
            </Alert>
          )}

          {hasReconstructions && itemScore < 0.5 && (
            <Stack direction="row" spacing={0.5} alignItems="center">
              <Tooltip
                title="Coverage is the share of the molecule RetroMol could match to known building blocks. Low or 0% coverage means most (or all) of the structure falls outside its rule set, e.g. an unusual scaffold or a modification the current rules don't recognize. A low coverage leads to primary sequences that may be sparse, incomplete, or empty."
                placement="bottom-start"
                arrow
              >
                <Stack direction="row" spacing={0.5} alignItems="center" sx={{ cursor: "help", width: "fit-content" }}>
                  <InfoOutlinedIcon fontSize="small" color="warning" />
                  <Typography variant="caption" color="text.secondary">
                    Low coverage: parts of the structure below may be sparse or empty
                  </Typography>
                </Stack>
              </Tooltip>
            </Stack>
          )}

          {dbMatchingSequences.length > 0 && (
            <Stack spacing={1}>
              <PrimarySequenceOverview
                data={dbMatchingSequences}
                selectedTags={selectedTags}
                onToggleMotif={handleToggleMotif}
              />
            </Stack>
          )}

          <Box sx={{ display: "flex", justifyContent: "flex-end" }}>
            <ExportImageButton
              targetRef={diagramRef}
              filename={`retromol-${(item.name || item.id).replace(/[^a-z0-9]+/gi, "-")}`}
              label="Download the diagram below as a PNG"
            />
          </Box>
          <Box
            ref={diagramRef}
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
                  smiles={hasReconstructions ? (data?.[0]?.tagged_input_smiles ?? "") : item.smiles}
                  size={300}
                  highlightAtoms={hasReconstructions ? highlightAtoms : []}
                />
              </ErrorBoundary>
            </Box>
            {hasReconstructions && (
              <>
                <AnnotatedArrow
                  annotation={
                    isUnordered
                      ? 'Parsed out motifs'
                      : allHaveBackbone
                      ? 'Linearization'
                      : 'Linearization and sequencing'
                  }
                />
                {allHaveBackbone && (
                  <>
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
                                smiles={reconstruction.tagged_backbone_smiles ?? ""}
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
                  </>
                )}
                <Box
                  sx={{
                    // Ordered: sized to its natural content width (not squeezed by
                    // the other flex:1 siblings) -- the full sequence should be
                    // readable without an extra inner scrollbar; the outer row above
                    // already scrolls horizontally for anything that doesn't fit.
                    // Unordered: given a bounded width instead, so the "bag" of
                    // motifs actually wraps onto multiple lines rather than
                    // stretching into one long unwrapped row.
                    flex: isUnordered ? '1 1 0' : '0 0 auto',
                    minWidth: 0,
                    maxWidth: isUnordered ? 480 : undefined,
                    display: 'flex',
                    justifyContent: 'center',
                    alignItems: 'center',
                    pr: 2
                  }}
                >
                  <Box
                    sx={{
                      flex: 1,
                      minWidth: 0,
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
                        minWidth: 0,
                        width: "100%",
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
              </>
            )}
          </Box>
          <DrawingAttribution library="smiles-drawer" sx={{ textAlign: "center" }} />
          {hasReconstructions && (
            <DescriptionBox
              title={'Explanation'}
              description={
                isUnordered
                  ? 'The input SMILES structure (left) is processed by RetroMol into non-overlapping building blocks (right), but no single biosynthetic order could be determined between them. They\'re listed in no particular order, not as a sequence. You can still highlight individual motifs by clicking them.'
                  : allHaveBackbone
                  ? 'The input SMILES structure (left) is processed by RetroMol into non-overlapping building blocks, and from these building blocks a linear backbone is reconstructed (middle). The primary sequence, a text representation of the linear backbone, is seen on the right. You can highlight individual motifs by clicking them in the primary sequence above. If a sequence was parsed wrong or incompletely, use "Edit sequences" to fix it by hand. Saved edits are kept with this compound and reused when you query it from the Discovery tab; "Revert" on a row discards the edit and restores the algorithm’s own parse.'
                  : 'The input SMILES structure (left) is processed by RetroMol into non-overlapping building blocks, and directly linearized into the primary sequence on the right (no reconstructed backbone structure is shown). You can highlight individual motifs by clicking them in the primary sequence above. If a sequence was parsed wrong or incompletely, use "Edit sequences" to fix it by hand. Saved edits are kept with this compound and reused when you query it from the Discovery tab; "Revert" on a row discards the edit and restores the algorithm’s own parse.'
              }
            />
          )}
        </Box>
      )}
    </DialogWindow>
  );
};