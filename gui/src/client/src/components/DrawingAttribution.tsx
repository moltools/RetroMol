import React from "react";
import Typography from "@mui/material/Typography";
import type { SxProps, Theme } from "@mui/material/styles";

// smiles-drawer's package.json is blocked from import by its own "exports" field, so
// this is kept in sync by hand with the resolved version in package-lock.json --
// scripts/checkSmilesDrawerVersion.js fails CI/lint if the two drift apart.
const SMILES_DRAWER_VERSION = "2.4.1";

type DrawingAttributionProps =
  | { library: "smiles-drawer"; sx?: SxProps<Theme> }
  // RDKit drawings are produced server-side (see gui/src/server/routes/rules.py), so
  // there's no local constant to hardcode -- the version is only known once the
  // caller has actually fetched a drawing and the server reported what rendered it.
  | { library: "rdkit"; version: string; sx?: SxProps<Theme> };

// A single caption attributing one or more structure/reaction drawings above it to the
// library that rendered them. Place once per diagram -- if a diagram contains multiple
// drawings from the same library (e.g. a compound plus its reconstructions), attribute
// the whole group once rather than repeating this per drawing.
export const DrawingAttribution: React.FC<DrawingAttributionProps> = (props) => {
  const label = props.library === "smiles-drawer" ? `SmilesDrawer v${SMILES_DRAWER_VERSION}` : `RDKit v${props.version}`;
  return (
    <Typography
      variant="caption"
      color="text.secondary"
      sx={{ display: "block", fontStyle: "italic", ...props.sx }}
    >
      Drawn with {label}
    </Typography>
  );
};
