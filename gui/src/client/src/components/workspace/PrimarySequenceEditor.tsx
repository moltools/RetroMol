import React from "react";
import Box from "@mui/material/Box";
import Chip from "@mui/material/Chip";
import CircularProgress from "@mui/material/CircularProgress";
import IconButton from "@mui/material/IconButton";
import Stack from "@mui/material/Stack";
import Tooltip from "@mui/material/Tooltip";
import Typography from "@mui/material/Typography";
import RestartAltIcon from "@mui/icons-material/RestartAlt";
import { Session, SessionItem } from "../../features/session/types";
import type { PrimarySequenceItem, Reconstruction } from "../../features/reconstruction/types";
import { saveEditedPrimarySequences, revertEditedPrimarySequence } from "../../features/reconstruction/api";
import { useNotifications } from "../NotificationProvider";
import { MotifName } from "../MotifName";
import { horizontalScrollSx } from "../../theme/scrollbarSx";
import { SequenceEditor, type SequenceBlock } from "./SequenceEditor";

export function blocksFromSequence(sequence: PrimarySequenceItem[]): SequenceBlock[] {
  return sequence.map(([name, tags]) => ({ id: crypto.randomUUID(), name, tags }));
}

export function sequenceFromBlocks(blocks: SequenceBlock[]): PrimarySequenceItem[] {
  return blocks.map((b) => [b.name, b.tags ?? []]);
}

// A missing override means "use whatever the algorithm just parsed" -- see the
// comment on CompoundItemSchema.editedPrimarySequences for why we never store a
// copy of the original.
export function savedSequenceFor(item: SessionItem, idx: number, original: PrimarySequenceItem[]): PrimarySequenceItem[] {
  if (item.kind !== "compound") return original;
  return item.editedPrimarySequences?.[String(idx)] ?? original;
}

// Shared editing state machine for a compound's primary sequence(s) -- used by both
// the full "View item" dialog and the inline expandable row in the Upload list, so
// edit/save/revert behave identically (and stay in sync) no matter where they're
// triggered from.
//
// `resetSignal` lets a caller force a fresh draft (e.g. the dialog re-seeds on every
// open, even for the same item); pass a stable value (or omit it) to only reset when
// `data`/`item.id` actually change, which is what lets in-progress edits survive
// unrelated re-renders (like another item's status ticking over).
export function usePrimarySequenceEditor(
  session: Session,
  setSession: React.Dispatch<React.SetStateAction<Session | null>>,
  item: SessionItem,
  data: Reconstruction[] | null | undefined,
  resetSignal?: unknown,
) {
  const { pushNotification } = useNotifications();
  const [editing, setEditing] = React.useState(false);
  const [drafts, setDrafts] = React.useState<Record<number, SequenceBlock[]>>({});
  const [saving, setSaving] = React.useState(false);
  const [revertingIdx, setRevertingIdx] = React.useState<number | null>(null);

  React.useEffect(() => {
    if (!data) return;
    const initial: Record<number, SequenceBlock[]> = {};
    data.forEach((reconstruction, idx) => {
      initial[idx] = blocksFromSequence(savedSequenceFor(item, idx, reconstruction.primary_sequence));
    });
    setDrafts(initial);
    setEditing(false);
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [data, item.id, resetSignal]);

  const isRowDirty = (idx: number, original: PrimarySequenceItem[]): boolean => {
    const draft = drafts[idx];
    if (!draft) return false;
    return JSON.stringify(sequenceFromBlocks(draft)) !== JSON.stringify(savedSequenceFor(item, idx, original));
  };

  const anyDirty = (data ?? []).some((reconstruction, idx) => isRowDirty(idx, reconstruction.primary_sequence));

  const handleCancelEdit = () => {
    if (data) {
      const initial: Record<number, SequenceBlock[]> = {};
      data.forEach((reconstruction, idx) => {
        initial[idx] = blocksFromSequence(savedSequenceFor(item, idx, reconstruction.primary_sequence));
      });
      setDrafts(initial);
    }
    setEditing(false);
  };

  const handleSaveAll = async () => {
    if (item.kind !== "compound" || !data) return;

    const overrides: Record<number, PrimarySequenceItem[]> = {};
    data.forEach((reconstruction, idx) => {
      if (isRowDirty(idx, reconstruction.primary_sequence)) {
        overrides[idx] = sequenceFromBlocks(drafts[idx] ?? []);
      }
    });
    if (Object.keys(overrides).length === 0) return;

    setSaving(true);
    try {
      const nextSession = await saveEditedPrimarySequences(session, item.id, overrides);
      setSession(() => nextSession);
      pushNotification("Saved edited primary sequence(s).", "success");
      setEditing(false);
    } catch (err) {
      const msg = err instanceof Error ? err.message : String(err);
      pushNotification(`Failed to save edited sequence(s): ${msg}`, "error");
    } finally {
      setSaving(false);
    }
  };

  const handleRevertRow = async (idx: number, original: PrimarySequenceItem[]) => {
    setDrafts((prev) => ({ ...prev, [idx]: blocksFromSequence(original) }));

    if (item.kind !== "compound" || !item.editedPrimarySequences?.[String(idx)]) return;

    setRevertingIdx(idx);
    try {
      const nextSession = await revertEditedPrimarySequence(session, item.id, idx);
      setSession(() => nextSession);
      pushNotification("Reverted to the algorithm-parsed sequence.", "success");
    } catch (err) {
      const msg = err instanceof Error ? err.message : String(err);
      pushNotification(`Failed to revert sequence: ${msg}`, "error");
    } finally {
      setRevertingIdx(null);
    }
  };

  return {
    editing,
    setEditing,
    drafts,
    setDrafts,
    saving,
    revertingIdx,
    anyDirty,
    isRowDirty,
    handleCancelEdit,
    handleSaveAll,
    handleRevertRow,
  };
}

export type PrimarySequenceEditorState = ReturnType<typeof usePrimarySequenceEditor>;

function PrimarySequenceChips({
  sequence,
  selectedTags,
  onToggleMotif,
}: {
  sequence: PrimarySequenceItem[];
  selectedTags: number[];
  onToggleMotif: (tags: number[]) => void;
}) {
  return (
    <Box sx={{ display: "flex", flexWrap: "nowrap", gap: 1, alignItems: "center", ...horizontalScrollSx }}>
      {sequence.map(([name, tags], idx) => {
        const isSelected = tags.length > 0 && tags.every((tag) => selectedTags.includes(tag));

        return (
          <Box
            key={`${name}-${idx}`}
            onClick={() => onToggleMotif(tags)}
            sx={{
              px: 1.25,
              py: 0.75,
              borderRadius: 1,
              border: "1px solid",
              borderColor: isSelected ? "primary.main" : "divider",
              bgcolor: isSelected ? "primary.main" : "background.paper",
              color: isSelected ? "primary.contrastText" : "text.primary",
              fontSize: "0.875rem",
              fontWeight: 500,
              cursor: "pointer",
              userSelect: "none",
              whiteSpace: "nowrap",
              flexShrink: 0,
              "&:hover": {
                borderColor: "primary.main",
                bgcolor: isSelected ? "primary.dark" : "action.hover",
              },
            }}
          >
            <MotifName name={name} />
          </Box>
        );
      })}
    </Box>
  );
}

// Renders one row per reconstruction: a label (+ "Edited" chip), and either the
// read-only chip row or the live SequenceEditor + per-row revert control.
export function PrimarySequenceRows({
  item,
  data,
  state,
  selectedTags,
  onToggleMotif,
  labelWidth = 130,
}: {
  item: SessionItem;
  data: Reconstruction[];
  state: PrimarySequenceEditorState;
  selectedTags: number[];
  onToggleMotif: (tags: number[]) => void;
  labelWidth?: number;
}) {
  const { editing, drafts, setDrafts, revertingIdx, isRowDirty, handleRevertRow } = state;

  return (
    <>
      {data.map((reconstruction, idx) => {
        const override = item.kind === "compound" ? item.editedPrimarySequences?.[String(idx)] : undefined;
        const draftBlocks = drafts[idx] ?? blocksFromSequence(override ?? reconstruction.primary_sequence);
        const dirty = isRowDirty(idx, reconstruction.primary_sequence);

        return (
          <Box
            key={`primary-sequence-row-${idx}`}
            sx={{ display: "flex", flexDirection: "row", alignItems: "center", gap: 1, width: "100%", minWidth: 0 }}
          >
            <Stack direction="column" spacing={0.25} sx={{ minWidth: labelWidth, flexShrink: 0 }}>
              <Typography variant="caption" sx={{ color: "text.secondary", fontWeight: 600 }}>
                {`primary sequence ${idx + 1}`}
              </Typography>
              {override && !editing && (
                <Chip
                  label="Edited"
                  size="small"
                  color="info"
                  variant="outlined"
                  sx={{ height: 18, fontSize: "0.65rem", width: "fit-content" }}
                />
              )}
            </Stack>

            {editing ? (
              <>
                <Box sx={{ flex: 1, minWidth: 0 }}>
                  <SequenceEditor
                    blocks={draftBlocks}
                    onChange={(blocks) => setDrafts((prev) => ({ ...prev, [idx]: blocks }))}
                    showProvenance
                    selectedTags={selectedTags}
                    onBlockClick={onToggleMotif}
                  />
                </Box>
                <Tooltip title="Revert to the algorithm-parsed sequence" arrow>
                  <span>
                    <IconButton
                      size="small"
                      disabled={!override && !dirty}
                      onClick={() => handleRevertRow(idx, reconstruction.primary_sequence)}
                    >
                      {revertingIdx === idx ? <CircularProgress size={16} /> : <RestartAltIcon fontSize="small" />}
                    </IconButton>
                  </span>
                </Tooltip>
              </>
            ) : (
              <Box sx={{ flex: 1, minWidth: 0 }}>
                <PrimarySequenceChips
                  sequence={override ?? reconstruction.primary_sequence}
                  selectedTags={selectedTags}
                  onToggleMotif={onToggleMotif}
                />
              </Box>
            )}
          </Box>
        );
      })}
    </>
  );
}
