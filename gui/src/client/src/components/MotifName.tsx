import React from "react";
import Box from "@mui/material/Box";

// Renders a monomer name like "A2^R" or "C^E2", superscripting the trailing stereo marker.
export function MotifName({ name }: { name: string }) {
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
