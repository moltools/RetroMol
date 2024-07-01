import React from "react";
import { Container, CssBaseline, Box } from "@mui/material";

const Page = ({ children, ...props }) => {
    return (
        <Container {...props}>
            <CssBaseline />
            <Box>
                {children}
            </Box>
        </Container>
    );
};

export default Page;