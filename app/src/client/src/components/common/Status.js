import React from "react";
import { Box } from "@mui/material";

const Status = ({ statusName, status }) => (
    <Box
        variant="contained"
        sx={{
            mb: 1,
            backgroundColor: status ? "#28c840" : "#ff5f57",
            transition: "background-color 0.5s ease",
            color: status ? "#006200" : "#990000",
            textAlign: "center",
            padding: "8px",
            borderRadius: "4px",
            cursor: "default",
            userSelect: "none",
            mx: 1
        }}
    >
        {status ? `${statusName}: Online` : `${statusName}: Offline`}
    </Box>
);

export default Status;