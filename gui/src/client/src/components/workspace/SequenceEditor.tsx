import React from "react";
import Box from "@mui/material/Box";
import IconButton from "@mui/material/IconButton";
import TextField from "@mui/material/TextField";
import Tooltip from "@mui/material/Tooltip";
import Autocomplete from "@mui/material/Autocomplete";
import DragIndicatorIcon from "@mui/icons-material/DragIndicator";
import CloseIcon from "@mui/icons-material/Close";
import AddIcon from "@mui/icons-material/Add";
import { useQuery } from "@tanstack/react-query";
import {
  DndContext,
  closestCenter,
  PointerSensor,
  KeyboardSensor,
  useSensor,
  useSensors,
  type DragEndEvent,
} from "@dnd-kit/core";
import {
  SortableContext,
  horizontalListSortingStrategy,
  useSortable,
  arrayMove,
  sortableKeyboardCoordinates,
} from "@dnd-kit/sortable";
import { CSS } from "@dnd-kit/utilities";
import { MotifName } from "../MotifName";
import { horizontalScrollSx } from "../../theme/scrollbarSx";
import { searchMonomerNames } from "../../features/discovery/api";
import type { MonomerNameOption } from "../../features/discovery/types";

// tags carries the source atom indices this block was mined from (if any). A
// block with no tags (freshly added, or otherwise hand-edited) has no atoms to
// highlight -- see `showProvenance` below.
export type SequenceBlock = { id: string; name: string; tags?: number[] };

type SequenceEditorProps = {
  blocks: SequenceBlock[];
  onChange: (blocks: SequenceBlock[]) => void;
  disabled?: boolean;
  // When true, visually distinguishes blocks that are still linked to real
  // parsed atoms (solid border, clickable) from ones that aren't (dashed
  // border) -- e.g. when this editor sits next to the structure it was mined
  // from and edits shouldn't be mistaken for algorithm output.
  showProvenance?: boolean;
  selectedTags?: number[];
  onBlockClick?: (tags: number[]) => void;
};

function SortableBlock({
  block,
  disabled,
  showProvenance,
  selected,
  onDelete,
  onClick,
}: {
  block: SequenceBlock;
  disabled?: boolean;
  showProvenance?: boolean;
  selected?: boolean;
  onDelete: (id: string) => void;
  onClick?: () => void;
}) {
  const { attributes, listeners, setNodeRef, transform, transition, isDragging } = useSortable({
    id: block.id,
    disabled,
  });

  const linked = (block.tags?.length ?? 0) > 0;
  const clickable = showProvenance && linked && !!onClick;

  const content = (
    <Box
      ref={setNodeRef}
      onClick={clickable ? onClick : undefined}
      sx={{
        transform: CSS.Transform.toString(transform),
        transition: transition ?? undefined,
        opacity: isDragging ? 0.5 : 1,
        display: "flex",
        alignItems: "center",
        gap: 0.5,
        px: 1,
        py: 0.75,
        borderRadius: 1,
        border: "1px solid",
        borderStyle: showProvenance && !linked ? "dashed" : "solid",
        borderColor: selected ? "primary.main" : "divider",
        bgcolor: selected ? "primary.main" : "background.paper",
        color: selected ? "primary.contrastText" : "text.primary",
        fontSize: "0.875rem",
        fontWeight: 500,
        whiteSpace: "nowrap",
        flexShrink: 0,
        cursor: clickable ? "pointer" : undefined,
        userSelect: clickable ? "none" : undefined,
      }}
    >
      <Box
        {...attributes}
        {...listeners}
        sx={{
          display: "flex",
          alignItems: "center",
          cursor: disabled ? "default" : "grab",
          color: selected ? "inherit" : "text.secondary",
          touchAction: "none",
        }}
      >
        <DragIndicatorIcon fontSize="small" />
      </Box>
      <Box component="span">
        <MotifName name={block.name} />
      </Box>
      {!disabled && (
        <IconButton
          size="small"
          onClick={(e) => {
            e.stopPropagation();
            onDelete(block.id);
          }}
          sx={{ p: 0.25, color: selected ? "inherit" : undefined }}
        >
          <CloseIcon fontSize="inherit" />
        </IconButton>
      )}
    </Box>
  );

  if (!showProvenance) return content;

  return (
    <Tooltip
      title={linked ? "Click to highlight the source atoms" : "Not linked to the original structure (added or edited by hand)"}
      arrow
    >
      {content}
    </Tooltip>
  );
}

function AddBlockControl({ disabled, onAdd }: { disabled?: boolean; onAdd: (name: string) => void }) {
  const [open, setOpen] = React.useState(false);
  const [inputValue, setInputValue] = React.useState("");
  const [selected, setSelected] = React.useState<MonomerNameOption | null>(null);

  const searchQuery = useQuery({
    queryKey: ["discoveryMonomerNames", inputValue.trim()],
    queryFn: ({ signal }) => searchMonomerNames(inputValue.trim(), 10, signal),
    enabled: open && inputValue.trim().length > 0,
    staleTime: 60_000,
  });

  const options = searchQuery.data ?? [];

  const reset = () => {
    setOpen(false);
    setInputValue("");
    setSelected(null);
  };

  const handleConfirm = () => {
    if (!selected) return;
    onAdd(selected.name);
    reset();
  };

  if (!open) {
    return (
      <IconButton
        size="small"
        disabled={disabled}
        onClick={() => setOpen(true)}
        sx={{ border: "1px dashed", borderColor: "divider", borderRadius: 1, flexShrink: 0 }}
      >
        <AddIcon fontSize="small" />
      </IconButton>
    );
  }

  return (
    <Box sx={{ display: "flex", alignItems: "center", gap: 0.5, flexShrink: 0 }}>
      {/* freeSolo is intentionally off: block names must resolve to a real monomer identity */}
      <Autocomplete<MonomerNameOption, false, false, false>
        size="small"
        sx={{ width: 180 }}
        options={options}
        loading={searchQuery.isFetching}
        getOptionLabel={(option) => option.name}
        isOptionEqualToValue={(option, value) => option.name === value.name}
        inputValue={inputValue}
        onInputChange={(_, value, reason) => {
          setInputValue(value);
          // MUI fires this right after onChange too (reason "reset", syncing the
          // input text to the newly-picked option's label) -- clearing `selected`
          // unconditionally here would immediately undo that selection. Only clear
          // it when the user is actually typing.
          if (reason === "input") {
            setSelected(null);
          }
        }}
        value={selected}
        onChange={(_, value) => setSelected(value)}
        renderInput={(params) => (
          <TextField {...params} autoFocus size="small" hiddenLabel placeholder="Block name..." />
        )}
      />
      <IconButton size="small" disabled={!selected} onClick={handleConfirm}>
        <AddIcon fontSize="small" />
      </IconButton>
      <IconButton size="small" onClick={reset}>
        <CloseIcon fontSize="small" />
      </IconButton>
    </Box>
  );
}

export const SequenceEditor: React.FC<SequenceEditorProps> = ({
  blocks,
  onChange,
  disabled = false,
  showProvenance = false,
  selectedTags = [],
  onBlockClick,
}) => {
  const sensors = useSensors(
    useSensor(PointerSensor, { activationConstraint: { distance: 4 } }),
    useSensor(KeyboardSensor, { coordinateGetter: sortableKeyboardCoordinates })
  );

  const handleDragEnd = (event: DragEndEvent) => {
    const { active, over } = event;
    if (!over || active.id === over.id) return;

    const oldIndex = blocks.findIndex((b) => b.id === active.id);
    const newIndex = blocks.findIndex((b) => b.id === over.id);
    if (oldIndex === -1 || newIndex === -1) return;

    onChange(arrayMove(blocks, oldIndex, newIndex));
  };

  const handleDelete = (id: string) => {
    onChange(blocks.filter((b) => b.id !== id));
  };

  // Requirement: adding always appends to the end. Manually added blocks carry
  // no tags -- there's no parsed atom range for them to link to.
  const handleAdd = (name: string) => {
    onChange([...blocks, { id: crypto.randomUUID(), name, tags: [] }]);
  };

  return (
    <DndContext sensors={sensors} collisionDetection={closestCenter} onDragEnd={handleDragEnd}>
      <Box sx={{ display: "flex", flexWrap: "nowrap", gap: 1, alignItems: "center", ...horizontalScrollSx }}>
        <SortableContext items={blocks.map((b) => b.id)} strategy={horizontalListSortingStrategy}>
          {blocks.map((block) => {
            const linked = (block.tags?.length ?? 0) > 0;
            const selected = linked && block.tags!.every((tag) => selectedTags.includes(tag));
            return (
              <SortableBlock
                key={block.id}
                block={block}
                disabled={disabled}
                showProvenance={showProvenance}
                selected={selected}
                onDelete={handleDelete}
                onClick={block.tags ? () => onBlockClick?.(block.tags!) : undefined}
              />
            );
          })}
        </SortableContext>
        {!disabled && <AddBlockControl disabled={disabled} onAdd={handleAdd} />}
      </Box>
    </DndContext>
  );
};
