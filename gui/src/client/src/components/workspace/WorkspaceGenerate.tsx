import React from "react";
import Box from "@mui/material/Box";
import Card from "@mui/material/Card";
import CardContent from "@mui/material/CardContent";
import Stack from "@mui/material/Stack";
import Typography from "@mui/material/Typography";
import Alert from "@mui/material/Alert";
import Button from "@mui/material/Button";
import CircularProgress from "@mui/material/CircularProgress";
import ContentCopyIcon from "@mui/icons-material/ContentCopy";
import { generateBackbone } from "../../features/reconstruction/api";
import type { GeneratedBackbone } from "../../features/reconstruction/types";
import { useNotifications } from "../NotificationProvider";
import { ErrorBoundary } from "../ErrorBoundary";
import { ExportImageButton } from "../ExportImageButton";
import SmilesDrawerContainer from "../SmilesDrawerContainer.js";
import { DrawingAttribution } from "../DrawingAttribution";
import { SequenceEditor, type SequenceBlock } from "./SequenceEditor";
import { PrimarySequenceChips } from "./PrimarySequenceEditor";

type HighlightAtom = [number, string];

export const WorkspaceGenerate: React.FC = () => {
  const { pushNotification } = useNotifications();

  const [blocks, setBlocks] = React.useState<SequenceBlock[]>([]);
  const [generating, setGenerating] = React.useState(false);
  const [result, setResult] = React.useState<GeneratedBackbone | null>(null);
  const [error, setError] = React.useState<string | null>(null);
  const [selectedTags, setSelectedTags] = React.useState<number[]>([]);
  const diagramRef = React.useRef<HTMLDivElement>(null);

  const canGenerate = blocks.length > 0 && !generating;

  const handleGenerate = async () => {
    if (!canGenerate) return;
    setGenerating(true);
    setError(null);
    setSelectedTags([]);

    try {
      const reconstruction = await generateBackbone(blocks.map((b) => b.name));
      setResult(reconstruction);
    } catch (err) {
      const msg = err instanceof Error ? err.message : String(err);
      setResult(null);
      setError(msg);
      pushNotification(`Failed to generate backbone: ${msg}`, "error");
    } finally {
      setGenerating(false);
    }
  };

  const handleToggleMotif = (tags: number[]) => {
    setSelectedTags((prev) => {
      const allSelected = tags.every((tag) => prev.includes(tag));
      if (allSelected) return prev.filter((tag) => !tags.includes(tag));
      return Array.from(new Set([...prev, ...tags]));
    });
  };

  const handleCopySmiles = async () => {
    if (!result?.backboneSmiles) return;
    try {
      await navigator.clipboard.writeText(result.backboneSmiles);
      pushNotification("Copied SMILES to clipboard", "success");
    } catch {
      pushNotification("Failed to copy SMILES", "error");
    }
  };

  // Memoized so it's a stable reference across re-renders that don't touch
  // selectedTags -- otherwise SmilesDrawerContainer redraws on every unrelated
  // parent re-render (same reasoning as DialogViewItem's highlightAtoms).
  const highlightAtoms = React.useMemo<HighlightAtom[]>(
    () => selectedTags.map((tag) => [tag, "#027bf3"]),
    [selectedTags]
  );

  const orientationTags = React.useMemo(() => {
    if (!result || result.primary_sequence.length < 2) return undefined;
    const [, startTags] = result.primary_sequence[0];
    const [, endTags] = result.primary_sequence[result.primary_sequence.length - 1];
    if (!startTags.length || !endTags.length) return undefined;
    return { startTags, endTags };
  }, [result]);

  return (
    <Box sx={{ width: "100%", mx: "auto", display: "flex", flexDirection: "column", gap: "16px" }}>
      <Card variant="outlined">
        <CardContent>
          <Typography component="h1" variant="subtitle1">
            Build a primary sequence
          </Typography>
          <Typography variant="body2" color="text.secondary" sx={{ mb: 1.5 }}>
            Add building blocks in biosynthetic order, then generate the linear backbone
            structure RetroMol's fusion chemistry would assemble from them.
          </Typography>

          {blocks.length === 0 ? (
            <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
              No blocks yet.
            </Typography>
          ) : null}

          <SequenceEditor blocks={blocks} onChange={setBlocks} disabled={generating} />

          <Stack direction="row" spacing={1.5} alignItems="center" sx={{ mt: 2 }}>
            <Button variant="contained" disabled={!canGenerate} onClick={handleGenerate}>
              Generate
            </Button>
            {generating && <CircularProgress size={20} />}
          </Stack>
        </CardContent>
      </Card>

      {(result || error) && (
        <Card variant="outlined">
          <CardContent>
            <Typography component="h1" variant="subtitle1">
              Generated backbone
            </Typography>

            {error && <Alert severity="error" sx={{ mt: 1 }}>{error}</Alert>}

            {result && (
              <Box sx={{ mt: 1.5 }}>
                {result.backbone_warning && (
                  <Alert severity="warning" sx={{ mb: 1.5 }}>
                    {result.backbone_warning}
                  </Alert>
                )}

                {result.tagged_backbone_smiles ? (
                  <>
                    <Box sx={{ display: "flex", justifyContent: "flex-end" }}>
                      <Stack direction="column" spacing={0} alignItems="flex-start">
                        <ExportImageButton
                          targetRef={diagramRef}
                          filename={`retromol-generated-${blocks.map((b) => b.name).join("-").replace(/[^a-z0-9]+/gi, "-")}`}
                          label="Download the diagram and sequence below as a PNG"
                        />
                        {result.backboneSmiles && (
                          <Button
                            size="small"
                            variant="text"
                            startIcon={<ContentCopyIcon fontSize="small" />}
                            onClick={handleCopySmiles}
                          >
                            Copy SMILES
                          </Button>
                        )}
                      </Stack>
                    </Box>
                    <Box ref={diagramRef} sx={{ display: "flex", flexDirection: "column", alignItems: "center", gap: 0.5, mb: 1.5 }}>
                      <ErrorBoundary
                        what="molecule structure"
                        fallback={
                          <Typography variant="caption" color="text.secondary">
                            Could not render this structure.
                          </Typography>
                        }
                      >
                        <SmilesDrawerContainer
                          identifier="generated-backbone"
                          smiles={result.tagged_backbone_smiles}
                          size={300}
                          highlightAtoms={highlightAtoms}
                          orientationTags={orientationTags}
                        />
                      </ErrorBoundary>
                      <DrawingAttribution library="smiles-drawer" />
                      <Box sx={{ display: "flex", justifyContent: "center", width: "100%", mt: 2 }}>
                        <PrimarySequenceChips
                          sequence={result.primary_sequence}
                          selectedTags={selectedTags}
                          onToggleMotif={handleToggleMotif}
                        />
                      </Box>
                    </Box>
                  </>
                ) : (
                  <Box sx={{ display: "flex", justifyContent: "center" }}>
                    <PrimarySequenceChips sequence={result.primary_sequence} selectedTags={[]} />
                  </Box>
                )}
              </Box>
            )}
          </CardContent>
        </Card>
      )}
    </Box>
  );
};
