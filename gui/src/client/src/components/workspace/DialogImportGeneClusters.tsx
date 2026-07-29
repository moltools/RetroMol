import React from "react";
import Stack from "@mui/material/Stack";
import Typography from "@mui/material/Typography";
import Button from "@mui/material/Button";
import MuiLink from "@mui/material/Link";
import Slider from "@mui/material/Slider";
import { DialogWindow } from "../DialogWindow";

// Matches ParasModel's own default (retromol_antismash.inference.model_paras).
const DEFAULT_PARAS_THRESHOLD = 0.1;

type DialogImportGeneClusterProps = {
  open: boolean;
  onClose: () => void;
  onImport: (files: File[], parasThreshold: number) => void;
}

export const DialogImportGeneCluster: React.FC<DialogImportGeneClusterProps> = ({
  open,
  onClose,
  onImport,
}) => {
  const [gbkFiles, setGbkFiles] = React.useState<File[]>([]);
  const [parasThreshold, setParasThreshold] = React.useState<number>(DEFAULT_PARAS_THRESHOLD);
  const canImport = gbkFiles.length > 0;

  const reset = () => {
    setGbkFiles([]);
    setParasThreshold(DEFAULT_PARAS_THRESHOLD);
  };

  const handleImport = () => {
    onImport(gbkFiles, parasThreshold);
    reset();
    onClose();
  }

  return (
    <DialogWindow
      open={open}
      onClose={onClose}
      title="Import gene clusters"
      dividers
      actions={[
        { label: "Cancel", variant: "text", color: "inherit", onClick: onClose },
        { label: "Clear", variant: "contained", color: "secondary", onClick: reset },
        { label: "Import", variant: "contained", color: "primary", onClick: handleImport, disabled: !canImport, autoFocus: true },
      ]}
    >
      <Stack spacing={2}>
        <Typography>
          Select one or more GenBank files (.gbk, .gb, .genbank) containing gene cluster data to import into your workspace. Make sure the files are&nbsp;
          <MuiLink href="https://antismash.secondarymetabolites.org/#!/start" target="_blank" rel="noopener noreferrer">
            antiSMASH
          </MuiLink>
          &nbsp;output files for best compatibility.
        </Typography>
        <Button variant="outlined" component="label">
          Choose files
          <input
            type="file"
            hidden
            multiple
            accept=".gb,.gbk,.genbank,application/genbank"
            onChange={(e) => setGbkFiles(Array.from(e.target.files || []))}
          />
        </Button>
        {gbkFiles.length > 0 && (
          <Typography variant="body2">
            {gbkFiles.length} file(s) selected
          </Typography>
        )}

        <Stack spacing={0.5}>
          <Typography variant="body2">
            PARAS substrate confidence threshold: {parasThreshold.toFixed(2)}
          </Typography>
          <Typography variant="caption" color="text.secondary">
            Minimum prediction probability required to call a substrate for an NRPS module. Lower surfaces more (lower-confidence) predictions; higher keeps only the most confident ones. Doesn't affect PKS modules, which are classified directly from domain anatomy.
          </Typography>
          <Slider
            value={parasThreshold}
            onChange={(_, value) => setParasThreshold(value as number)}
            min={0}
            max={1}
            step={0.01}
            marks={[
              { value: 0, label: "0" },
              { value: DEFAULT_PARAS_THRESHOLD, label: "default" },
              { value: 1, label: "1" },
            ]}
            valueLabelDisplay="auto"
            valueLabelFormat={(value) => value.toFixed(2)}
            sx={{ maxWidth: 360 }}
          />
        </Stack>
      </Stack>
    </DialogWindow>
  )
}