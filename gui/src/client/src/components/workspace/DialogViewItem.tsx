import React from "react";
import Alert from "@mui/material/Alert";
import CircularProgress from "@mui/material/CircularProgress";
import Box from "@mui/material/Box";
import Typography from "@mui/material/Typography";
import { SessionItem } from "../../features/session/types";
import { DialogWindow } from "../DialogWindow";
import SmilesDrawerContainer from "../SmilesDrawerContainer.js";

export type PrimarySequenceItem = [string, number[]];

export type Reconstruction = {
  tagged_input_smiles: string;
  tagged_backbone_smiles: string;
  primary_sequence: PrimarySequenceItem[];
};

type HighlightAtom = [number, string];

type DialogViewItemProps = {
  sessionId: string;
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

function PrimarySequence({
  sequence,
  selectedTags,
  onToggleMotif,
}: {
  sequence: PrimarySequenceItem[];
  selectedTags: number[];
  onToggleMotif: (tags: number[]) => void;
}) {
  return (
    <Box
      sx={{
        display: "flex",
        flexWrap: "nowrap",
        gap: 1,
        justifyContent: "center",
        alignItems: "center",
      }}
    >
      {sequence.map(([name, tags], idx) => {
        const isSelected =
          tags.length > 0 && tags.every((tag) => selectedTags.includes(tag));

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
              "&:hover": {
                borderColor: "primary.main",
                bgcolor: isSelected ? "primary.dark" : "action.hover",
              },
            }}
          >
            {name}
          </Box>
        )
      })}
    </Box>
  )
}

export const DialogViewItem: React.FC<DialogViewItemProps> = ({
  sessionId,
  item,
  open,
  onClose,
}) => {
  const isCompound = item.kind === "compound"; // there are only two types: "compound" and "cluster"

  const [loading, setLoading] = React.useState<boolean>(false);
  const [error, setError] = React.useState<string | null>(null);
  const [data, setData] = React.useState<Reconstruction | null>(null);
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

  const highlightAtoms: HighlightAtom[] = selectedTags.map((tag) => [tag, "#027bf3"])

  // Clear selection when dialog/item changes
  React.useEffect(() => {
    if (open) {
      setSelectedTags([]);
    }
  }, [open, item.id]);

  React.useEffect(() => {
    if (!open) { return; };
    if (!isCompound) { return; };
    const fetchData = async () => {
      setLoading(true);
      setError(null);
      try {
        const response = await fetch("/api/reconstructCompound", {
          method: "POST",
          headers: {
            "Content-Type": "application/json",
          },
          body: JSON.stringify({
            sessionId,
            itemId: item.id,
          }),
        });
        if (!response.ok) { throw new Error(`Error fetching data for linear view: ${response.statusText}`); };
        const data = await response.json();
        setData(data.data);
      } catch (err: any) {
        setError(err.message || "Unknown error");
      } finally {
        setLoading(false);
      }
    }
    fetchData();
  }, [open, isCompound, sessionId, item.id]);

  return (
    <DialogWindow
      open={open}
      onClose={onClose}
      title="View item"
      dividers
      actions={[
        { label: "Close", variant: "text", color: "inherit", onClick: onClose },
      ]}
      maxWidth={"lg"}
    >
      {!isCompound && (
        <Alert severity="info" sx={{ mb: 2 }}>
          Viewing is only available for compounds.
        </Alert>
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
              <SmilesDrawerContainer
                identifier={`smiles-drawer-${sessionId}-${item.id}-full`}
                smiles={data?.tagged_input_smiles ?? ""}
                size={300}
                highlightAtoms={highlightAtoms}
              />
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
              <SmilesDrawerContainer
                identifier={`smiles-drawer-${sessionId}-${item.id}-preprocessed`}
                smiles={data?.tagged_backbone_smiles ?? ""}
                size={300}
                highlightAtoms={highlightAtoms}
              />
            </Box>
            <AnnotatedArrow annotation={'Sequencing'} />
            <Box
              sx={{
                flex: 1,
                display: 'flex',
                justifyContent: 'center',
                alignItems: 'center',
                pr: 2
              }}
            >
              <PrimarySequence
                sequence={data?.primary_sequence ?? []}
                selectedTags={selectedTags}
                onToggleMotif={handleToggleMotif}
              />
            </Box>
          </Box>
          <DescriptionBox
            title={'Explanation'}
            description={'The input SMILES (right) is processed by RetroMol into non-overlapping building blocks, and from these building blocks a linear backbone is reconstructed (middle). The primary sequence, a text representation of the linear backbone, is seen on the right. You can highlight individual motifs by clicking them in the primary sequence above.'}
          />
        </Box>
      )}
    </DialogWindow>
  );
};