import React from "react";
import Box from "@mui/material/Box";
import Chip from "@mui/material/Chip";
import CircularProgress from "@mui/material/CircularProgress";
import Tooltip from "@mui/material/Tooltip";
import Typography from "@mui/material/Typography";
import RestartAltIcon from "@mui/icons-material/RestartAlt";
import { Session, SessionItem } from "../../features/session/types";
import type { PrimarySequenceItem, Reconstruction } from "../../features/reconstruction/types";
import { saveEditedPrimarySequences, revertEditedPrimarySequence } from "../../features/reconstruction/api";
import { useNotifications } from "../NotificationProvider";
import { MotifHoverCard } from "../MotifHoverCard";
import { MotifName } from "../MotifName";
import { horizontalScrollSx } from "../../theme/scrollbarSx";
import { SequenceEditor, type SequenceBlock } from "./SequenceEditor";
import { MinimalIconButton } from "../MinimalIconButton";

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
  ordered = true,
}: {
  sequence: PrimarySequenceItem[];
  selectedTags: number[];
  // Omit this to render a plain, non-interactive chip row -- e.g. the inline
  // preview in the Upload list, which has no molecule view alongside it for a
  // highlight to make sense against. Only the "View item" dialog wires this up.
  onToggleMotif?: (tags: number[]) => void;
  // False renders this as an unordered bag of motifs instead of a sequence: no
  // connecting line, no fixed reading order (wraps instead of scrolling as one
  // row) -- see Reconstruction.ordered. Must stay visually distinct from the
  // ordered case so it's never mistaken for a real primary sequence.
  ordered?: boolean;
}) {
  const selectable = !!onToggleMotif;

  return (
    <Box sx={ordered ? { ...horizontalScrollSx } : undefined}>
      <Box
        sx={{
          display: "flex",
          flexWrap: ordered ? "nowrap" : "wrap",
          gap: 0.5,
          alignItems: "center",
          position: "relative",
          width: ordered ? "max-content" : "100%",

          // The connecting "string" -- only for an actual sequence; an unordered
          // bag of motifs has no line, that's the whole point.
          ...(ordered && {
            "&::before": {
              content: '""',
              position: "absolute",
              left: 0,
              right: 0,
              top: "50%",
              height: "2px",
              backgroundColor: "divider",
              zIndex: 0,
            },
          }),
        }}
      >
        {sequence.map(([name, tags], idx) => {
          const linked = tags.length > 0;
          const isSelected = selectable && linked && tags.every((tag) => selectedTags.includes(tag));

          const chip = (
            <Box
              key={`${name}-${idx}`}
              onClick={selectable ? () => onToggleMotif!(tags) : undefined}
              sx={{
                px: 1.25,
                py: 0.75,
                borderRadius: 1,
                border: "1px solid",
                borderStyle: linked ? "solid" : "dashed",
                borderColor: isSelected ? "primary.main" : "divider",
                bgcolor: isSelected ? "primary.main" : "background.paper",
                color: isSelected ? "primary.contrastText" : "text.primary",
                fontSize: "0.875rem",
                fontWeight: 500,
                cursor: selectable ? "pointer" : "default",
                userSelect: "none",
                whiteSpace: "nowrap",
                flexShrink: 0,
                position: "relative",
                zIndex: 1,
                ...(selectable && {
                  "&:hover": {
                    borderColor: "primary.main",
                  },
                }),
              }}
            >
              <MotifName name={name} />
            </Box>
          );

          // Same rule as the drag-and-drop editor's blocks: only claim it's
          // clickable when it actually is (see SortableBlock in SequenceEditor).
          const hint = !linked
            ? "Not linked to the original structure (added or edited by hand)"
            : selectable
            ? "Click to highlight the source atoms"
            : undefined;

          return (
            <MotifHoverCard key={`${name}-${idx}-hover`} name={name} hint={hint}>
              {chip}
            </MotifHoverCard>
          );
        })}
      </Box>
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
  labelWidth,
}: {
  item: SessionItem;
  data: Reconstruction[];
  state: PrimarySequenceEditorState;
  selectedTags: number[];
  // Omit to disable click-to-highlight entirely (see PrimarySequenceChips) --
  // used for both the read-only chip row and the drag-and-drop editor's blocks.
  onToggleMotif?: (tags: number[]) => void;
  labelWidth?: number;
}) {
  const { editing, drafts, setDrafts, revertingIdx, isRowDirty, handleRevertRow } = state;

  const labelColumn = labelWidth ? `${labelWidth}px` : "max-content";

  return (
    <>
      {data.map((reconstruction, idx) => {
        const override =
          item.kind === "compound"
            ? item.editedPrimarySequences?.[String(idx)]
            : undefined;

        const draftBlocks =
          drafts[idx] ??
          blocksFromSequence(override ?? reconstruction.primary_sequence);

        const dirty = isRowDirty(idx, reconstruction.primary_sequence);
        const ordered = reconstruction.ordered !== false;

        return (
          <Box
            key={`primary-sequence-row-${idx}`}
            sx={{
              display: "grid",
              gridTemplateColumns: `${labelColumn} minmax(0, 1fr)`,
              alignItems: "center",
              columnGap: 1,
              width: "100%",
              minWidth: 0,
            }}
          >
            {/* Label column */}
            <Box
              sx={{
                display: "flex",
                alignItems: "center",
                gap: 0.75,
                minWidth: 0,
                overflow: "hidden",
                transform: ordered ? "translateY(-5px)" : "none",
              }}
            >
              <Typography
                variant="caption"
                sx={{
                  color: "text.secondary",
                  fontWeight: 600,
                  whiteSpace: "nowrap",
                  overflow: "hidden",
                  textOverflow: "ellipsis",
                  minWidth: 0,
                }}
              >
                {ordered ? `primary sequence ${idx + 1}` : "parsed motifs (unordered)"}
              </Typography>

              {override && !editing && (
                <Chip
                  label="Edited"
                  size="small"
                  color="info"
                  variant="outlined"
                  sx={{
                    height: 18,
                    fontSize: "0.65rem",
                    flexShrink: 0,
                  }}
                />
              )}
            </Box>

            {/* Sequence column */}
            <Box sx={{ minWidth: 0 }}>
              {editing && ordered ? (
                <Box sx={{ display: "flex", alignItems: "center", gap: 1 }}>
                  <Box sx={{ flex: 1, minWidth: 0 }}>
                    <SequenceEditor
                      blocks={draftBlocks}
                      onChange={(blocks) =>
                        setDrafts((prev) => ({
                          ...prev,
                          [idx]: blocks,
                        }))
                      }
                      showProvenance
                      selectedTags={selectedTags}
                      onBlockClick={onToggleMotif}
                    />
                  </Box>

                  <Tooltip title="Revert to the algorithm-parsed sequence" arrow>
                    <span>
                      <MinimalIconButton
                        size="small"
                        disabled={!override && !dirty}
                        onClick={() =>
                          handleRevertRow(idx, reconstruction.primary_sequence)
                        }
                        sx={{ transform: "translateY(-6px)"}}
                      >
                        {revertingIdx === idx ? (
                          <CircularProgress size={16} />
                        ) : (
                          <RestartAltIcon fontSize="small" />
                        )}
                      </MinimalIconButton>
                    </span>
                  </Tooltip>
                </Box>
              ) : (
                <PrimarySequenceChips
                  sequence={override ?? reconstruction.primary_sequence}
                  selectedTags={selectedTags}
                  onToggleMotif={onToggleMotif}
                  ordered={ordered}
                />
              )}
            </Box>
          </Box>
        );
      })}
    </>
  );
}
