import React from "react";
import Stack from "@mui/material/Stack";
import Typography from "@mui/material/Typography";
import Button from "@mui/material/Button";
import MuiLink from "@mui/material/Link";
import { DialogWindow } from "../DialogWindow";

type DialogImportGeneClusterProps = {
  open: boolean;
  onClose: () => void;
  onImport: (files: File[]) => void;
}

export const DialogImportGeneCluster: React.FC<DialogImportGeneClusterProps> = ({
  open,
  onClose,
  onImport,
}) => {
  const [gbkFiles, setGbkFiles] = React.useState<File[]>([]);
  const canImport = gbkFiles.length > 0;

  const reset = () => setGbkFiles([]);

  const handleImport = () => {
    onImport(gbkFiles);
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
      </Stack>
    </DialogWindow>
  )
}