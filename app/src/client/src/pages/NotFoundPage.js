import React from "react";
import { Link as RouterLink } from "react-router-dom";
import { Typography, Button, Container } from "@mui/material";

const NotFoundPage = () => {
    return (
        <Container maxWidth="sm" style={{ textAlign: 'center', marginTop: '50px' }}>
            <Typography variant="h2" gutterBottom>
                404 - Page Not Found
            </Typography>
            <Typography variant="body1" gutterBottom>
                The page you are looking for might have been removed, had its name changed, or is temporarily unavailable.
            </Typography>
            <Button variant="contained" color="primary" component={RouterLink} to="/" sx={{ marginTop: '1rem' }}>
                Go Back to Landing Page
            </Button>
        </Container>
    );
};

export default NotFoundPage;