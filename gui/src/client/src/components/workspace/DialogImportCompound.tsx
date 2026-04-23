import React from "react";
import Chip from "@mui/material/Chip";
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

type CompoundOption = {
  name: string;
  smiles: string;
  databaseName: string;
  databaseIdentifier: string;
}

type DialogImportCompoundProps = {
  open: boolean;
  onClose: () => void;
  onImportSingle: (compound: { name: string; smiles: string; matchStereochemistry: boolean }) => void;
  onImportBatch: (file: File, matchStereochemistry: boolean) => void;
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

  // Stereochemistry toggle
  const [matchStereochemistry, setMatchStereochemistry] = React.useState<boolean>(true);

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
  };

  const handleImport = () => {
    if (mode === "single") {
      onImportSingle({
        name: compoundName.trim(),
        smiles: compoundSmiles.trim(),
        matchStereochemistry,
      });
    } else if (batchFile) {
      onImportBatch(batchFile, matchStereochemistry);
    }
    reset();
    onClose();
  };

  async function searchCompoundByName(q: string) {
    const params = new URLSearchParams({q, limit: "10"});
    const res = await fetch(`/api/searchCompound?${params.toString()}`);
    if (!res.ok) { throw new Error(`Search failed: ${res.status}`); };
    return await res.json();
  };

  // Debounced search when user types a compound name
  React.useEffect(() => {
    if (!open) return;
    const q = compoundName.trim();
    if (!q) {
      setOptions([]);
      return;
    };

    const handle = setTimeout(async () => {
      setLoading(true);
      try {
        const res = await searchCompoundByName(q);
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
        Begin typing a compound name to see autocomplete suggestions from our database. Selecting one will automatically fill in a valid SMILES.
      </Typography>

      <Stack direction="column" spacing={0.5}>
        <FormControlLabel
          control={
            <Switch
              checked={matchStereochemistry}
              onChange={(e) => setMatchStereochemistry(e.target.checked)}
            />
          }
          label={matchStereochemistry ? "Match stereochemistry" : "Ignore stereochemistry"}
        />

        <FormControlLabel
          control={
            <Switch
              checked={mode === "batch"}
              onChange={(e) => setMode(e.target.checked ? "batch" : "single")}
            />
          }
          label={mode === "batch" ? "Batch import from file" : "Single compound import"}
        />
      </Stack>

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
            renderOption={(props, option) => {
              const { key, ...optionProps } = props as React.HTMLAttributes<HTMLLIElement> & { key: React.Key };

              if (typeof option === "string") {
                return (
                  <li key={key} {...optionProps}>
                    {option}
                  </li>
                );
              };

              return (
                <li key={key} {...optionProps}>
                  <Stack direction="row" spacing={1} alignItems="center">
                    <Chip
                      label={option.databaseName}
                      size="small"
                      variant="outlined"
                    />
                    <Typography variant="body2">
                      {option.name}
                    </Typography>
                  </Stack>
                </li>
              )
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
  );
};