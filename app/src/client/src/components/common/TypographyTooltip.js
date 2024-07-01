import React from "react";
import { Box, Tooltip } from "@mui/material";
import { Help as HelpIcon } from '@mui/icons-material';

const TypographyTooltip = ({ title, fontSize, children }) => {
    return (
        <Box sx={{ display: "flex", gap: 1 }}>
            {children}
            <Tooltip title={title} arrow>
                <HelpIcon sx={{ alignSelf: "center", fontSize: fontSize || 20, marginBottom: "6px" }} />
            </Tooltip>
        </Box>
    );
};

export default TypographyTooltip;