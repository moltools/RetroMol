import React from "react";
import Box from "@mui/material/Box";
import packageJson from "../../package.json";

export const BuildVersion: React.FC = () => {
  return (
    <Box
      sx={{
        display: "flex",
        flexDirection: "column",
        justifyContent: "center",
        alignItems: "center",
        fontSize: "0.8rem",
        color: "#888",
        m: 1,
      }}
    >
      <span>
        build: {packageJson.version}
      </span>
    </Box>
  )
}