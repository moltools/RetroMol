import React from "react";
import Box from "@mui/material/Box";
import { isLinkToken } from "../features/reconstruction/types";

// Stereo markers mxn.yml appends after "^": R/S (PK alpha/beta carbons), E/Z (PK
// alkene geometry), L/D (NRPS amino acid alpha carbons, e.g. "alanine^L").
const STEREO_MARKER = /^\^[A-Za-z]+$/;

// Renders a monomer name like "A2^R" or "C^E2" or "alanine^L", superscripting the
// trailing stereo marker(s). The link token (joining two merged primary sequence
// paths, see LINK_TOKEN) isn't a real building block -- render it as a muted glyph
// instead of a name, everywhere a monomer name would otherwise be rendered (sequence
// chips, alignment grid cells, ...).
export function MotifName({ name }: { name: string }) {
  if (isLinkToken(name)) {
    return (
      <Box component="span" sx={{ color: "text.disabled", fontWeight: 400 }} title="Joins two merged primary sequence paths">
        &#8942;
      </Box>
    );
  }

  const parts = name.split(/(\^[A-Za-z]+)/g);

  return (
    <>
      {parts.map((part, i) => {
        if (STEREO_MARKER.test(part)) {
          return (
            <Box
              key={i}
              component="sup"
              sx={{
                fontSize: "0.7em",
                lineHeight: 0,
                verticalAlign: "super",
              }}
            >
              {part.slice(1)}
            </Box>
          );
        }

        return <React.Fragment key={i}>{part}</React.Fragment>;
      })}
    </>
  );
}
