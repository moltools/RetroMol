import React from "react";
import Typography from "@mui/material/Typography";
import type { SxProps, Theme } from "@mui/material/styles";

// Neither version is importable at runtime -- smiles-drawer's package.json is blocked
// by its own "exports" field, and RDKit drawings are produced server-side (see
// gui/src/server/routes/rules.py). Both are kept in sync by hand with the resolved
// version in package-lock.json and pyproject.toml's `rdkit==` pin, respectively.
const RDKIT_VERSION = "2025.9.1";
const SMILES_DRAWER_VERSION = "2.4.1";

type DrawingLibrary = "smiles-drawer" | "rdkit";

const LIBRARY_LABEL: Record<DrawingLibrary, string> = {
  "smiles-drawer": `SmilesDrawer v${SMILES_DRAWER_VERSION}`,
  rdkit: `RDKit v${RDKIT_VERSION}`,
};

// A single caption attributing one or more structure/reaction drawings above it to the
// library that rendered them. Place once per diagram -- if a diagram contains multiple
// drawings from the same library (e.g. a compound plus its reconstructions), attribute
// the whole group once rather than repeating this per drawing.
export const DrawingAttribution: React.FC<{ library: DrawingLibrary; sx?: SxProps<Theme> }> = ({ library, sx }) => {
  return (
    <Typography
      variant="caption"
      color="text.secondary"
      sx={{ display: "block", fontStyle: "italic", ...sx }}
    >
      Drawn with {LIBRARY_LABEL[library]}
    </Typography>
  );
};
