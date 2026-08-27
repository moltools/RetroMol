import React from "react";
import Box from "@mui/material/Box";
import CircularProgress from "@mui/material/CircularProgress";
import Chip from "@mui/material/Chip";
import Stack from "@mui/material/Stack";
import Tooltip from "@mui/material/Tooltip";
import Typography from "@mui/material/Typography";
import RestartAltIcon from "@mui/icons-material/RestartAlt";
import { ClusterItem, Session, SessionItem } from "../../features/session/types";
import type { ClusterModule, ClusterPrimarySequence } from "../../features/clusters/types";
import { saveEditedClusterPrimarySequence, revertEditedClusterPrimarySequence } from "../../features/clusters/api";
import type { PrimarySequenceItem } from "../../features/reconstruction/types";
import { blocksFromSequence, sequenceFromBlocks, PrimarySequenceChips } from "./PrimarySequenceEditor";
import { SequenceEditor, type SequenceBlock } from "./SequenceEditor";
import { AnnotatedArrow } from "./AnnotatedArrow";
import { MinimalIconButton } from "../MinimalIconButton";
import { useNotifications } from "../NotificationProvider";
import { horizontalScrollSx } from "../../theme/scrollbarSx";

// Mirrors the PKS extender-unit classification in retromol_antismash/modules.py
// (PKSModule.substrate) for the tooltip only -- the chip's own label comes from
// the region's resolved primary sequence (see ClusterReadoutDiagram), which
// already accounts for the specific extender digit/stereo suffix a module may
// have resolved (pmp.yml's "pks.extender"/"pks.kr_stereochemistry" axes), not
// just the reduction letter these anatomy flags alone can tell you.
function pksReductionLetter(anatomy: Extract<ClusterModule, { type: "PKS" }>["anatomy"]): string {
  const KR = anatomy.has_active_KR;
  const DH = anatomy.has_active_DH;
  const ER = anatomy.has_active_ER;
  const effDH = DH && KR;
  const effER = ER && KR && DH;

  if (!KR) return "A";
  if (!effDH) return "B";
  if (!effER) return "C";
  return "D";
}

// A single gene/domain module box in the "genes & domains" row -- selectable the
// same way a SequenceEditor block or PrimarySequenceChips entry is, so clicking a
// readout token can highlight the module(s) it came from and vice versa.
function ModuleBox({
  module,
  selected,
  onClick,
}: {
  module: ClusterModule;
  selected: boolean;
  onClick?: () => void;
}) {
  const isNRPS = module.type === "NRPS";
  const substrateName = isNRPS ? module.predicted_substrate?.name ?? null : null;
  const hasConfidentSubstrate = isNRPS && !!substrateName && substrateName !== "unknown";

  const tooltipLines = [
    `${module.type} module ${module.module_index_in_gene + 1} in ${module.gene_id} (${module.gene_strand} strand)`,
    `domains: ${module.present_domains.join(", ") || "none"}`,
    !isNRPS ? `reduction level: PKS_${pksReductionLetter(module.anatomy)}` : null,
    isNRPS && module.predicted_substrate?.score != null
      ? `prediction confidence: ${(module.predicted_substrate.score * 100).toFixed(1)}%`
      : null,
    onClick ? "Click to highlight the matching readout token" : null,
  ].filter(Boolean).join("\n");

  return (
    <Tooltip title={<span style={{ whiteSpace: "pre-line" }}>{tooltipLines}</span>} arrow>
      <Box
        onClick={onClick}
        sx={{
          px: 1.25,
          py: 0.75,
          borderRadius: 1,
          border: "1px solid",
          borderColor: selected ? "primary.main" : (isNRPS ? (hasConfidentSubstrate ? "success.main" : "divider") : "info.main"),
          bgcolor: selected ? "primary.main" : (isNRPS ? (hasConfidentSubstrate ? "success.main" : "background.paper") : "info.main"),
          color: selected ? "primary.contrastText" : (isNRPS ? (hasConfidentSubstrate ? "success.contrastText" : "text.primary") : "white"),
          fontSize: "0.8rem",
          whiteSpace: "nowrap",
          flexShrink: 0,
          display: "flex",
          flexDirection: "column",
          alignItems: "center",
          lineHeight: 1.3,
          minWidth: 72,
          cursor: onClick ? "pointer" : "default",
          userSelect: "none",
          "&:hover": onClick ? { borderColor: "primary.main" } : undefined,
        }}
      >
        <Typography variant="caption" sx={{ opacity: 0.75, fontSize: "0.65rem" }}>
          {module.type} · {module.gene_id}
        </Typography>
        <Typography variant="body2" sx={{ fontWeight: 600 }}>
          {isNRPS ? (hasConfidentSubstrate ? substrateName : "unknown") : `PKS_${pksReductionLetter(module.anatomy)}`}
        </Typography>
      </Box>
    </Tooltip>
  );
}

// The "genes & domains" row -- one selectable box per module, in the same
// biosynthetic order as the readout row below it (see ClusterPrimarySequence's
// index alignment, features/clusters/types.ts).
function ModulesRow({
  modules,
  selectedIndices,
  onToggleModule,
}: {
  modules: ClusterModule[];
  selectedIndices: number[];
  onToggleModule?: (index: number) => void;
}) {
  return (
    <Box sx={{ display: "flex", flexWrap: "nowrap", gap: 1, alignItems: "center", ...horizontalScrollSx }}>
      {modules.map((module, idx) => (
        <ModuleBox
          key={`${module.gene_id}-${idx}`}
          module={module}
          selected={selectedIndices.includes(idx)}
          onClick={onToggleModule ? () => onToggleModule(idx) : undefined}
        />
      ))}
    </Box>
  );
}

// Shared editing state machine for a gene cluster's per-region primary sequence --
// mirrors PrimarySequenceEditor.tsx's usePrimarySequenceEditor, but keyed by
// region index against ClusterItem.editedPrimarySequences instead of
// CompoundItem.editedPrimarySequences.
export function useClusterPrimarySequenceEditor(
  session: Session,
  setSession: React.Dispatch<React.SetStateAction<Session | null>>,
  item: SessionItem,
  data: ClusterPrimarySequence[] | null | undefined,
  resetSignal?: unknown,
) {
  const { pushNotification } = useNotifications();
  const [editing, setEditing] = React.useState(false);
  const [drafts, setDrafts] = React.useState<Record<number, SequenceBlock[]>>({});
  const [saving, setSaving] = React.useState(false);
  const [revertingIdx, setRevertingIdx] = React.useState<number | null>(null);

  const savedSequenceFor = (idx: number, original: PrimarySequenceItem[]): PrimarySequenceItem[] => {
    if (item.kind !== "cluster") return original;
    return (item as ClusterItem).editedPrimarySequences?.[String(idx)] ?? original;
  };

  React.useEffect(() => {
    if (!data) return;
    const initial: Record<number, SequenceBlock[]> = {};
    data.forEach((region, idx) => {
      initial[idx] = blocksFromSequence(savedSequenceFor(idx, region.primary_sequence));
    });
    setDrafts(initial);
    setEditing(false);
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [data, item.id, resetSignal]);

  const isRowDirty = (idx: number, original: PrimarySequenceItem[]): boolean => {
    const draft = drafts[idx];
    if (!draft) return false;
    return JSON.stringify(sequenceFromBlocks(draft)) !== JSON.stringify(savedSequenceFor(idx, original));
  };

  const anyDirty = (data ?? []).some((region, idx) => isRowDirty(idx, region.primary_sequence));

  const handleCancelEdit = () => {
    if (data) {
      const initial: Record<number, SequenceBlock[]> = {};
      data.forEach((region, idx) => {
        initial[idx] = blocksFromSequence(savedSequenceFor(idx, region.primary_sequence));
      });
      setDrafts(initial);
    }
    setEditing(false);
  };

  const handleSaveAll = async () => {
    if (item.kind !== "cluster" || !data) return;

    const dirtyIndices = data
      .map((region, idx) => idx)
      .filter((idx) => isRowDirty(idx, data[idx].primary_sequence));
    if (dirtyIndices.length === 0) return;

    setSaving(true);
    try {
      let nextSession = session;
      for (const idx of dirtyIndices) {
        nextSession = await saveEditedClusterPrimarySequence(nextSession, item.id, idx, sequenceFromBlocks(drafts[idx] ?? []));
      }
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

    if (item.kind !== "cluster" || !(item as ClusterItem).editedPrimarySequences?.[String(idx)]) return;

    setRevertingIdx(idx);
    try {
      const nextSession = await revertEditedClusterPrimarySequence(session, item.id, idx);
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
    savedSequenceFor,
    handleCancelEdit,
    handleSaveAll,
    handleRevertRow,
  };
}

export type ClusterPrimarySequenceEditorState = ReturnType<typeof useClusterPrimarySequenceEditor>;

// One region: a "genes & domains" row, an arrow, and the readout (primary
// sequence) row -- either read-only chips or, while editing, a live
// SequenceEditor. Clicking a chip on either row highlights its counterpart on
// the other, via `selectedIndices`/`onToggleIndex` (shared with the sibling row).
function RegionDiagram({
  item,
  region,
  regionIndex,
  state,
  selectedIndices,
  onToggleIndex,
}: {
  item: SessionItem;
  region: ClusterPrimarySequence;
  regionIndex: number;
  state: ClusterPrimarySequenceEditorState;
  selectedIndices: number[];
  onToggleIndex: (index: number) => void;
}) {
  const { editing, drafts, setDrafts, revertingIdx, isRowDirty, savedSequenceFor, handleRevertRow } = state;

  const override = item.kind === "cluster" ? (item as ClusterItem).editedPrimarySequences?.[String(regionIndex)] : undefined;
  const draftBlocks = drafts[regionIndex] ?? blocksFromSequence(override ?? region.primary_sequence);
  const dirty = isRowDirty(regionIndex, region.primary_sequence);

  const onToggleMotif = (tags: number[]) => onToggleIndex(tags[0]);

  return (
    <Box>
      <Stack direction="row" spacing={0.75} alignItems="center" sx={{ mb: 0.5 }}>
        <Typography variant="caption" sx={{ color: "text.secondary", fontWeight: 600 }}>
          {region.id} ({region.modules.length} module{region.modules.length === 1 ? "" : "s"})
        </Typography>
        {override && !editing && (
          <Chip label="edited" size="small" color="info" variant="outlined" sx={{ fontSize: "0.65rem" }} />
        )}
      </Stack>

      {region.modules.length === 0 ? (
        <Typography variant="body2" color="text.secondary">
          No NRPS/PKS modules were detected in this region.
        </Typography>
      ) : (
        // flex: '0 0 auto' throughout (not the default shrinkable flex item) -- every
        // section keeps its natural content width, so none of their own internal
        // horizontal scrollbars (ModulesRow's, PrimarySequenceChips'/SequenceEditor's)
        // ever has to activate. This whole row overflows instead, and the single
        // outer diagram container (DialogViewItem.tsx's clusterDiagramRef) scrolls it
        // as one unit -- same pattern as a compound's structure/backbone/sequence
        // diagram.
        <Box sx={{ display: "flex", alignItems: "center", gap: 2, width: "max-content" }}>
          <Box sx={{ flex: "0 0 auto" }}>
            <ModulesRow modules={region.modules} selectedIndices={selectedIndices} onToggleModule={onToggleIndex} />
          </Box>

          {/* Extra breathing room left+right of the arrow, beyond the flex gap already
              separating it from the boxes on either side. */}
          <Box sx={{ px: 1.5, flex: "0 0 auto" }}>
            <AnnotatedArrow annotation="Readout" />
          </Box>

          <RegionReadoutControls
            override={override}
            draftBlocks={draftBlocks}
            editing={editing}
            dirty={dirty}
            regionIndex={regionIndex}
            region={region}
            revertingIdx={revertingIdx}
            setDrafts={setDrafts}
            handleRevertRow={handleRevertRow}
            selectedIndices={selectedIndices}
            onToggleMotif={onToggleMotif}
            fitContent
          />
        </Box>
      )}
    </Box>
  );
}

// The editor-or-chips + revert control shared by RegionDiagram (the full gene->readout
// view in the "View item" dialog) and ClusterPrimarySequenceRows (the compound-style
// primary-sequence-only row used in the Upload tab's inline expand).
function RegionReadoutControls({
  override,
  draftBlocks,
  editing,
  dirty,
  regionIndex,
  region,
  revertingIdx,
  setDrafts,
  handleRevertRow,
  selectedIndices,
  onToggleMotif,
  fitContent = false,
}: {
  override: PrimarySequenceItem[] | undefined;
  draftBlocks: SequenceBlock[];
  editing: boolean;
  dirty: boolean;
  regionIndex: number;
  region: ClusterPrimarySequence;
  revertingIdx: number | null;
  setDrafts: React.Dispatch<React.SetStateAction<Record<number, SequenceBlock[]>>>;
  handleRevertRow: (idx: number, original: PrimarySequenceItem[]) => void;
  selectedIndices: number[];
  onToggleMotif: (tags: number[]) => void;
  // True for RegionDiagram (the full gene->readout view): the readout keeps its
  // natural content width instead of shrinking to fit, so PrimarySequenceChips'/
  // SequenceEditor's own internal scrollbar never has to activate -- the row as a
  // whole (genes + arrow + readout) overflows instead, letting the single outer
  // diagram container (see DialogViewItem.tsx's clusterDiagramRef) be the one and
  // only scroll unit, same as a compound's structure/backbone/sequence diagram.
  // False (default) for ClusterPrimarySequenceRows (the Upload tab's inline row,
  // sequence-only, no genes/arrow to combine with) -- there, shrinking to fit and
  // scrolling internally is the correct, compound-inline-row-matching behavior.
  fitContent?: boolean;
}) {
  const sectionSx = fitContent ? { flex: "0 0 auto" as const } : { flex: 1, minWidth: 0 };

  return (
    <Box sx={{ ...sectionSx, display: "flex", alignItems: "center", gap: 1 }}>
      {editing ? (
        <>
          <Box sx={sectionSx}>
            <SequenceEditor
              blocks={draftBlocks}
              onChange={(blocks) => setDrafts((prev) => ({ ...prev, [regionIndex]: blocks }))}
              showProvenance
              selectedTags={selectedIndices}
              onBlockClick={onToggleMotif}
            />
          </Box>
          <Tooltip title="Revert to the algorithm-parsed sequence" arrow>
            <span>
              <MinimalIconButton
                size="small"
                disabled={!override && !dirty}
                onClick={() => handleRevertRow(regionIndex, region.primary_sequence)}
                sx={{ transform: "translateY(-6px)"}}
              >
                {revertingIdx === regionIndex ? <CircularProgress size={16} /> : <RestartAltIcon fontSize="small" />}
              </MinimalIconButton>
            </span>
          </Tooltip>
        </>
      ) : (
        // mb cancels out PrimarySequenceChips' own reserved scrollbar space (see
        // horizontalScrollSx's pb) so it doesn't throw off alignItems: "center"
        // against the label chip -- same fix as PrimarySequenceOverview's row
        // (PrimarySequenceEditor.tsx) for a compound's primary sequence.
        <Box sx={{ ...sectionSx, pr: 2, mb: -1.5 }}>
          <PrimarySequenceChips
            sequence={override ?? region.primary_sequence}
            selectedTags={selectedIndices}
            onToggleMotif={onToggleMotif}
          />
        </Box>
      )}
    </Box>
  );
}

// Full "genes/domains -> readout" diagram for every region in a gene cluster,
// with click-to-highlight linking a readout token to the module(s) it came
// from. Wrap this in a ref and pass that ref to ExportImageButton to let a user
// download it as a PNG (see DialogViewItem.tsx for the pattern used for a
// compound's own structure/backbone/sequence diagram).
export function ClusterReadoutDiagram({
  item,
  data,
  state,
}: {
  item: SessionItem;
  data: ClusterPrimarySequence[];
  state: ClusterPrimarySequenceEditorState;
}) {
  const [selectedIndices, setSelectedIndices] = React.useState<number[]>([]);

  // Selection is per-region in spirit (a module index only makes sense within its
  // own region), but since only one region's chips are ever clicked at a time in
  // practice, a single shared selection list keeps the state simple; toggling
  // clears cleanly per click regardless of which region it came from.
  const handleToggleIndex = (index: number) => {
    setSelectedIndices((prev) => (prev.includes(index) ? prev.filter((i) => i !== index) : [...prev, index]));
  };

  if (data.length === 0) {
    return (
      <Typography variant="body2" color="text.secondary">
        No antiSMASH regions were found in this file.
      </Typography>
    );
  }

  return (
    <Stack spacing={3}>
      {data.map((region, idx) => (
        <RegionDiagram
          key={`${region.id}-${idx}`}
          item={item}
          region={region}
          regionIndex={idx}
          state={state}
          selectedIndices={selectedIndices}
          onToggleIndex={handleToggleIndex}
        />
      ))}
    </Stack>
  );
}

// Compound-style primary-sequence-only view of a gene cluster's regions -- a
// "primary sequence" tag + the readout chip row (or, while editing, the
// SequenceEditor) per region, same shape as PrimarySequenceOverview/
// PrimarySequenceRows for a compound. No genes/domains track, no arrow -- that's
// ClusterReadoutDiagram's job in the full "View item" dialog; this is for the
// Upload tab's inline expandable row.
export function ClusterPrimarySequenceRows({
  item,
  data,
  state,
}: {
  item: SessionItem;
  data: ClusterPrimarySequence[];
  state: ClusterPrimarySequenceEditorState;
}) {
  const { editing, drafts, setDrafts, revertingIdx, isRowDirty, handleRevertRow } = state;

  if (data.length === 0) {
    return (
      <Typography variant="body2" color="text.secondary">
        No antiSMASH regions were found in this file.
      </Typography>
    );
  }

  return (
    <Stack spacing={1.5}>
      {data.map((region, idx) => {
        const override = item.kind === "cluster" ? (item as ClusterItem).editedPrimarySequences?.[String(idx)] : undefined;
        const draftBlocks = drafts[idx] ?? blocksFromSequence(override ?? region.primary_sequence);
        const dirty = isRowDirty(idx, region.primary_sequence);

        return (
          <Box key={`${region.id}-${idx}`} sx={{ display: "flex", alignItems: "center", gap: 1.5, minWidth: 0 }}>
            <Tooltip title={`${region.id} contains ${region.modules.length} biosynthetic module${region.modules.length === 1 ? "" : "s"}.`} arrow>
              <Chip label="primary sequence" size="small" color="info" variant="outlined" sx={{ fontSize: "0.7rem", flexShrink: 0, cursor: "help" }} />
            </Tooltip>
            {override && !editing && (
              <Chip label="edited" size="small" color="info" variant="outlined" sx={{ fontSize: "0.7rem", flexShrink: 0 }} />
            )}
            <Box sx={{ flex: 1, minWidth: 0 }}>
              <RegionReadoutControls
                override={override}
                draftBlocks={draftBlocks}
                editing={editing}
                dirty={dirty}
                regionIndex={idx}
                region={region}
                revertingIdx={revertingIdx}
                setDrafts={setDrafts}
                handleRevertRow={handleRevertRow}
                selectedIndices={[]}
                onToggleMotif={() => {}}
              />
            </Box>
          </Box>
        );
      })}
    </Stack>
  );
}
