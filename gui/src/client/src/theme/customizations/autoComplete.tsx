import { Theme, Components } from "@mui/material/styles";

export const autoCompleteCustomizations: Components<Theme> = {
  MuiAutocomplete: {
    styleOverrides: {
      // customize the dropdown container (paper) to limit its max height
      paper: ({ theme }) => ({
        maxHeight: 200,
        boxShadow: theme.shadows[2],
      }),
      // clearIndicator should never have a background or border
      clearIndicator: {
        padding: 0,
        width: 24,
        height: 24,
        background: "none",
        border: "none",
        "&:hover": {
          background: "none",
          border: "none",
        },
      },
    },
  },
}
