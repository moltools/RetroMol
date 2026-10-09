import type { SxProps, Theme } from "@mui/material/styles";

// A horizontal scrollbar that stays visible whenever there's overflow (rather
// than the OS's auto-hide-on-idle overlay scrollbar), with enough reserved
// space below the content that the bar never overlaps it.
export const horizontalScrollSx: SxProps<Theme> = {
  overflowX: "auto",
  overflowY: "hidden",
  pb: 1.5,
  scrollbarWidth: "thin",
  scrollbarColor: (theme) => `${theme.palette.action.disabled} transparent`,
  "&::-webkit-scrollbar": {
    height: 8,
  },
  "&::-webkit-scrollbar-track": {
    background: "transparent",
  },
  "&::-webkit-scrollbar-thumb": {
    backgroundColor: (theme) => theme.palette.action.disabled,
    borderRadius: 4,
  },
};
