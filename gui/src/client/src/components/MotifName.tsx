import React from "react";
import Box from "@mui/material/Box";
import { isLinkToken } from "../features/reconstruction/types";

// Renders a monomer name like "A2^R" or "C^E2", superscripting the trailing stereo marker.
// The link token (joining two merged primary sequence paths, see LINK_TOKEN) isn't a real
// building block -- render it as a muted glyph instead of a name, everywhere a monomer
// name would otherwise be rendered (sequence chips, alignment grid cells, ...).
export function MotifName({ name }: { name: string }) {
  if (isLinkToken(name)) {
    return (
      <Box component="span" sx={{ color: "text.disabled", fontWeight: 400 }} title="Joins two merged primary sequence paths">
        &#8942;
      </Box>
    );
  }

  const parts = name.split(/(\^[SREZ])/g);

  return (
    <>
      {parts.map((part, i) => {
        if (part === "^S" || part === "^R" || part === "^E" || part === "^Z") {
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
