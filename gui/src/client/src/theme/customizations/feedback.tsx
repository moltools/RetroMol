import { Theme, alpha, Components } from "@mui/material/styles";
import { gray, orange } from "../themePrimitives";

export const feedbackCustomizations: Components<Theme> = {
  MuiAlert: {
    styleOverrides: {
      root: ({ theme }) => ({
        borderRadius: 10,
        backgroundColor: orange[100],
        color: (theme.vars || theme).palette.text.primary,
        border: `1px solid ${alpha(orange[300], 0.5)}`,
        "& .MuiAlert-icon": {
          color: orange[500],
        },
        ...theme.applyStyles("dark", {
          backgroundColor: `${alpha(orange[900], 0.5)}`,
          border: `1px solid ${alpha(orange[800], 0.5)}`,
        }),
      }),
    },
  },
  MuiDialog: {
    styleOverrides: {
      root: ({ theme }) => ({
        "& .MuiDialog-paper": {
          borderRadius: "10px",
          border: "1px solid",
          borderColor: (theme.vars || theme).palette.divider,
          // Dialog content leans on subtle divider-colored lines/borders (e.g. the
          // primary-sequence connecting line) that barely register against the
          // app's near-black default paper background -- lighten just the dialog
          // surface in dark mode so those low-contrast details stay legible.
          ...theme.applyStyles("dark", {
            backgroundColor: "hsl(220, 20%, 13%)",
          }),
        },
      }),
    },
  },
  MuiLinearProgress: {
    styleOverrides: {
      root: ({ theme }) => ({
        height: 8,
        borderRadius: 8,
        backgroundColor: gray[200],
        ...theme.applyStyles("dark", {
          backgroundColor: gray[800],
        }),
      }),
    },
  },
}
