import React from "react";
import Box from "@mui/material/Box";
import Button from "@mui/material/Button";
import FormControlLabel from "@mui/material/FormControlLabel";
import Stack from "@mui/material/Stack";
import Switch from "@mui/material/Switch";
import TextField from "@mui/material/TextField";
import Typography from "@mui/material/Typography";
import Autocomplete from "@mui/material/Autocomplete";
import { useNotifications } from "../NotificationProvider";
import { DialogWindow } from "../DialogWindow";
import { runQuery } from "../../features/query/api";

type CompoundOption = {
  name: string;
  smiles: string;
}

type DialogImportCompoundProps = {
  open: boolean;
  onClose: () => void;
  onImportSingle: (compound: { name: string; smiles: string }) => void;
  onImportBatch: (file: File) => void;
}

export const DialogImportCompound: React.FC<DialogImportCompoundProps> = ({
  open,
  onClose,
  onImportSingle,
  onImportBatch,
}) => {
  const { pushNotification } = useNotifications();

  const [mode, setMode] = React.useState<"single" | "batch">("single");
  const [compoundName, setCompoundName] = React.useState<string>("");
  const [compoundSmiles, setCompoundSmiles] = React.useState<string>("");
  const [batchFile, setBatchFile] = React.useState<File | null>(null);

  // Autocomplete state
  const [options, setOptions] = React.useState<CompoundOption[]>([]);
  const [loading, setLoading] = React.useState<boolean>(false);

  const canImport =
    (mode === "single" &&
      compoundSmiles.trim().length > 0 &&
      compoundName.trim().length > 0) ||
    (mode === "batch" && batchFile !== null)

  const reset = () => {
    setMode("single");
    setCompoundName("");
    setCompoundSmiles("");
    setBatchFile(null);
    setOptions([]);
  }

  const handleImport = () => {
    if (mode === "single") {
      onImportSingle({
        name: compoundName.trim(),
        smiles: compoundSmiles.trim(),
      });
    } else if (batchFile) {
      onImportBatch(batchFile);
    }
    reset();
    onClose();
  }

  // Debounced search when user types a compound name
  React.useEffect(() => {
    if (!open) return;
    const q = compoundName.trim();
    if (!q) {
      setOptions([]);
      return;
    }

    const handle = setTimeout(async () => {
      setLoading(true);
      try {
        const res = await runQuery({
          name: "search_compound_by_name",
          params: { q },
          paging: { limit: 10 }, // size of suggestion list
        });

        const rows = (res.rows || []) as CompoundOption[];
        setOptions(rows);
      } catch (err) {
        pushNotification("Error searching compounds: " + (err as Error).message, "error");
        setOptions([]);
      } finally {
        setLoading(false);
      }
    }, 250); // 250ms debounce

    return () => clearTimeout(handle);
  }, [compoundName, open]);

  return (
    <DialogWindow
      open={open}
      onClose={onClose}
      title="Import compounds"
      dividers
      actions={[
        { label: "Cancel", variant: "text", color: "inherit", onClick: () => onClose() },
        { label: "Clear", variant: "contained", color: "secondary", onClick: reset },
        { label: "Import", variant: "contained", color: "primary", onClick: handleImport, disabled: !canImport, autoFocus: true },
      ]}
    >

    <Stack spacing={2}>
      <Typography>
        Enter a single compound identifier & SMILES, or upload a CSV/TSV that contains a column called "name" and a column called "smiles".
      </Typography>

      <FormControlLabel
        control={
          <Switch
            checked={mode === "batch"}
            onChange={(e) => setMode(e.target.checked ? "batch" : "single")}
          />
        }
        label={mode === "batch" ? "Batch import from file" : "Single compound import"}
      />

      {/* Single compound input */}
      {mode === "single" && (
        <Stack spacing={2}>
          <Autocomplete<CompoundOption, false, false, true>
            freeSolo
            loading={loading}
            options={options}
            getOptionLabel={(option) => typeof option === "string" ? option : option.name}
            // What goes in the input box
            inputValue={compoundName}
            onInputChange={(_, value) => {
              setCompoundName(value);
              if (!value) {
                setCompoundSmiles("");
              }
            }}
            // What happens when user selects an option
            onChange={(_, value) => {
              if (!value) {
                setCompoundName("");
                setCompoundSmiles("");
              } else if (typeof value === "string") {
                setCompoundName(value);
                // no SMILES from DB in this case
              } else {
                setCompoundName(value.name);
                setCompoundSmiles(value.smiles); // auto-fill from DB
              }
            }}
            renderInput={(params) => (
              <TextField
                {...params}
                placeholder="Start typing compound name..."
                size="small"
                hiddenLabel
                fullWidth
              />
            )}
          />

          {/* SMILES is still editable */}
          <TextField
            value={compoundSmiles}
            onChange={(e) => setCompoundSmiles(e.target.value)}
            placeholder="SMILES"
            size="small"
            hiddenLabel
            fullWidth
          />
        </Stack>
      )}

      {/* Batch compound file input */}
      {mode === "batch" && (
        <Stack spacing={2}>
          <Button variant="outlined" component="label">
            Choose CSV/TSV file
            <input
              type="file"
              hidden
              accept=".csv,.tsv,text/csv,text/tab-separated-values,text/plain"
              onChange={(e) => setBatchFile(e.target.files?.[0] || null)}
            />
          </Button>

          {batchFile && (
            <Box sx={{ fontSize: "0.875rem", color: "text.secondary" }}>
              Selected file: {batchFile.name}
              <br />
              <em>Expected columns: name, smiles</em>
            </Box>
          )}
        </Stack>
      )}

      </Stack>
    </DialogWindow>
  )
}