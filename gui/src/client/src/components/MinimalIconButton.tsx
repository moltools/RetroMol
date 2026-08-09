import React from "react";
import IconButton, { IconButtonProps } from "@mui/material/IconButton";

// A borderless, paddingless IconButton for inline "icon next to text" affordances
// (e.g. a rename pencil sitting right next to a label) where a normal IconButton's
// padding/hover halo would look out of place. Forwards its ref so it still works
// as the direct child of a Tooltip.
export const MinimalIconButton = React.forwardRef<HTMLButtonElement, IconButtonProps>(
  ({ sx, size = "small", ...props }, ref) => (
    <IconButton
      ref={ref}
      size={size}
      sx={[
        {
          p: 0,
          m: 0,
          minWidth: 0,
          minHeight: 0,
          width: "auto",
          height: "auto",
          border: "none !important",
          borderRadius: 0,
          outline: "none",
          background: "transparent !important",
          boxShadow: "none !important",

          "&:hover": {
            background: "transparent !important",
          },

          "&:focus": {
            background: "transparent !important",
            outline: "none",
          },
        },
        ...(Array.isArray(sx) ? sx : [sx]),
      ]}
      {...props}
    />
  )
);

MinimalIconButton.displayName = "MinimalIconButton";
