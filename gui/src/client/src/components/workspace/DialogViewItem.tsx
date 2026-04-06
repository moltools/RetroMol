import React from "react";
import Alert from "@mui/material/Alert";
import Stack from "@mui/material/Stack";
import Typography from "@mui/material/Typography";
import FormControl from "@mui/material/FormControl";
import Select from "@mui/material/Select";
import MenuItem from "@mui/material/MenuItem";
import MuiLink from "@mui/material/Link";
import CircularProgress from "@mui/material/CircularProgress";
import { useNotifications } from "../NotificationProvider";
import { DialogWindow } from "../DialogWindow";
import { SessionItem } from "../../features/session/types";
import { SvgViewer } from "../SvgViewer";
import { drawCompoundItem, drawGeneClusterItem } from "../../features/drawing/api";

type DialogViewItemProps = {
  open: boolean;
  item?: SessionItem | null;
  onClose: () => void;
}

export const DialogViewItem: React.FC<DialogViewItemProps> = ({
  open,
  item,
  onClose,
}) => {
  const { pushNotification } = useNotifications();
  const [loading, setLoading] = React.useState<boolean>(false);
  const [errorMsg, setErrorMsg] = React.useState<string | null>(null);

  const [initializedItemId, setInitializedItemId] = React.useState<string | null>(null);
  const [selectPrimarySequenceId, setSelectPrimarySequenceId] = React.useState<string>("");
  const [svg, setSvg] = React.useState<string | null>(null);

  const generateSvg = React.useCallback(async (primarySequenceId: string) => {
    if (!item) {
      setSvg(null);
      setErrorMsg(null);
      return;
    }

    // Get Item primarySequence for primarySequenceId
    const primarySequence = item.primarySequences?.find(seq => seq.id === primarySequenceId);
    if (!primarySequence) {
      pushNotification("Selected primary sequence not found in item.", "error");
      setSvg(null);
      setErrorMsg(null);
      return;
    }

    setLoading(true);

    try {
      if (item.kind === "compound") {
        const taggedParentSmiles = item.taggedSmiles;

        // Check if taggedParentSmiles is available
        if (!taggedParentSmiles) {
          pushNotification("No tagged SMILES available for this compound item.", "error");
          setSvg(null);
          return;
        }

        // Call drawing API
        const drawingSvg = await drawCompoundItem(
          taggedParentSmiles,
          primarySequence
        );
        setSvg(drawingSvg);
        setErrorMsg(null);
      } else if (item.kind === "gene_cluster") {
        // Call drawing API
        const drawingSvg = await drawGeneClusterItem(item.fileContent || "");
        setSvg(drawingSvg);
        setErrorMsg(null);
      } else {
        const errorMsg = "SVG drawing not supported for this item type.";
        pushNotification(errorMsg, "error");
        setSvg(null);
        setErrorMsg(errorMsg);
        return;
      }
    } catch (error) {
      let errorMsg = "Error generating SVG drawing";
      const errorBody = (error as any)?.body as string | undefined
      if (errorBody) {
        try {
          const parsed = JSON.parse(errorBody);
          if (typeof parsed?.error === "string") errorMsg = `${errorMsg}: ${parsed.error}`;
        } catch (e) {
          errorMsg = `${errorMsg}.`
        }
      }
      pushNotification(errorMsg, "error");
      setSvg(null);
      setErrorMsg(errorMsg);
    } finally {
      setLoading(false);
    }
  }, [item]);

  const initializePrimarySequence = React.useCallback(() => {
    if (item && item.primarySequences && item.primarySequences.length > 0) {
      const firstSeqId = item.primarySequences[0].id;
      setSelectPrimarySequenceId(firstSeqId);
      generateSvg(firstSeqId);
    } else {
      setSelectPrimarySequenceId("");
      setSvg(null);
      setErrorMsg(null);
    }
  }, [item, generateSvg]);

  // Initialize primary sequence selection when item changes
  // Also avoid re-initializing if the same item is passed again
  React.useEffect(() => {
    if (!item) {
      setInitializedItemId(null);
      setSelectPrimarySequenceId("");
      setSvg(null);
      setErrorMsg(null);
      return;
    }

    // No redraw if same item
    if (initializedItemId === item.id) {
      return;
    }

    // New item is passed
    setInitializedItemId(item.id);
    initializePrimarySequence();
  }, [item, initializedItemId, initializePrimarySequence]);

  const handlePrimarySequenceChange = (sequenceId: string) => {
    setSelectPrimarySequenceId(sequenceId);
    generateSvg(sequenceId);
  }

  return (
    <DialogWindow
      open={open}
      onClose={onClose}
      title="View item"
      dividers
      maxWidth="xl"
      actions={[
        { label: "Close", variant: "text", color: "inherit", onClick: onClose },
      ]}
    >
      { loading ? (
        <Stack direction="row" justifyContent="center" alignItems="center" height={400}>
          <CircularProgress />
        </Stack>
      ) : errorMsg ? (
        <Alert severity="error">{errorMsg}</Alert>
      ) : !item ? (
        <Typography variant="body1">
          No item selected.
        </Typography>
      ) : svg === null || svg.length === 0 ? (
        <Typography variant="body1">
          No SVG drawing available for this item.
        </Typography>
      ) : item.kind === "compound" ? (
        <Stack direction="column" gap={1} alignItems="flex-start">
          <Typography variant="body1">
            To get started, select a primary sequence to map onto the input compound structure.
            A downloadable SVG will be generated showing the mapping below the selector upon successful mapping.
            All structures in this view are drawn using&nbsp;
            <MuiLink href="https://research.wur.nl/en/publications/pikachu-a-python-based-informatics-kit-for-analysing-chemical-uni/" target="_blank" rel="noopener noreferrer">
              PIKAChU
            </MuiLink>
            &nbsp;. Please cite PIKAChU if you use these drawings in your work.
          </Typography>
          <FormControl fullWidth size="small">
            <Select
              label="Primary sequence"
              value={selectPrimarySequenceId}
              onChange={(e) => handlePrimarySequenceChange(e.target.value)}
              disabled={!item.primarySequences || item.primarySequences.length === 0}
            >
              {item.primarySequences && item.primarySequences.map((seq) => (
                <MenuItem key={seq.id} value={seq.id}>
                  {seq.name || seq.id}
                </MenuItem>
              ))}
            </Select>
          </FormControl>
          {svg && (
            <SvgViewer
              svg={svg}
              onZoomChange={() => {}}
              onElementClick={() => {}}
              height={600}
            />
          )}
        </Stack>
      ) : item.kind === "gene_cluster" ? (
        <Stack direction="column" spacing={2} alignItems="flex-start">
          <Typography variant="body1">
            <MuiLink href="https://research.wur.nl/en/publications/raichu-automating-the-visualisation-of-natural-product-biosynthes/" target="_blank" rel="noopener noreferrer">
              RAIChU
            </MuiLink>
            &nbsp;is used to generate the gene cluster visualization below.
            The SVG rendered below is included for informative purposes and does not reflect the exact encoding mechanism of the gene clustering by RetroMol.
            This viewer serves as a wrapper around the RAIChU SVG generation API. For more information on RAIChU, please refer to the link provided.
            Substrate predictions for non-ribosomal peptide (NRP) A-domains and polyketide synthase (PKS) acyltransferase (AT)-domains are taken directly from the&nbsp;
            <MuiLink href="https://antismash.secondarymetabolites.org/#!/start" target="_blank" rel="noopener noreferrer">
              antiSMASH
            </MuiLink>
            &nbsp;output.
            This wrapper viewer around RAIChU currently does not use any of the PARAS substrate specificity predictions provided and used by RetroMol for similarity searches.
            Additionally, this wrapper view around RAIChU only provided a visualization of the full region readout, not of individual candidate clusters found within the region.
            Please cite RAIChU if you use these drawings in your work.
          </Typography>
          {svg && (
            <SvgViewer
              svg={svg}
              onZoomChange={() => {}}
              onElementClick={() => {}}
              height={600}
            />
          )}
        </Stack>
      ) : (
        <Typography variant="body1">
          No preview available for this item type.
        </Typography>
      )}

    </DialogWindow>
  )
}