import React from "react";
import Box from "@mui/material/Box";
import IconButton from "@mui/material/IconButton";
import TextField from "@mui/material/TextField";
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

export type SequenceBlock = { id: string; name: string };

type SequenceEditorProps = {
  blocks: SequenceBlock[];
  onChange: (blocks: SequenceBlock[]) => void;
  disabled?: boolean;
};

function SortableBlock({
  block,
  disabled,
  onDelete,
}: {
  block: SequenceBlock;
  disabled?: boolean;
  onDelete: (id: string) => void;
}) {
  const { attributes, listeners, setNodeRef, transform, transition, isDragging } = useSortable({
    id: block.id,
    disabled,
  });

  return (
    <Box
      ref={setNodeRef}
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
        borderColor: "divider",
        bgcolor: "background.paper",
        fontSize: "0.875rem",
        fontWeight: 500,
        whiteSpace: "nowrap",
        flexShrink: 0,
      }}
    >
      <Box
        {...attributes}
        {...listeners}
        sx={{
          display: "flex",
          alignItems: "center",
          cursor: disabled ? "default" : "grab",
          color: "text.secondary",
          touchAction: "none",
        }}
      >
        <DragIndicatorIcon fontSize="small" />
      </Box>
      <Box component="span">
        <MotifName name={block.name} />
      </Box>
      {!disabled && (
        <IconButton size="small" onClick={() => onDelete(block.id)} sx={{ p: 0.25 }}>
          <CloseIcon fontSize="inherit" />
        </IconButton>
      )}
    </Box>
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
        onInputChange={(_, value) => {
          setInputValue(value);
          setSelected(null);
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

export const SequenceEditor: React.FC<SequenceEditorProps> = ({ blocks, onChange, disabled = false }) => {
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

  // Requirement: adding always appends to the end.
  const handleAdd = (name: string) => {
    onChange([...blocks, { id: crypto.randomUUID(), name }]);
  };

  return (
    <DndContext sensors={sensors} collisionDetection={closestCenter} onDragEnd={handleDragEnd}>
      <Box sx={{ display: "flex", flexWrap: "nowrap", gap: 1, alignItems: "center", ...horizontalScrollSx }}>
        <SortableContext items={blocks.map((b) => b.id)} strategy={horizontalListSortingStrategy}>
          {blocks.map((block) => (
            <SortableBlock key={block.id} block={block} disabled={disabled} onDelete={handleDelete} />
          ))}
        </SortableContext>
        {!disabled && <AddBlockControl disabled={disabled} onAdd={handleAdd} />}
      </Box>
    </DndContext>
  );
};
