import React from 'react';
import { Link } from 'react-router-dom';
import { Box, Typography, Button } from '@mui/material';

/**
 * Component to display a 404 error page.
 * 
 * @returns {React.ReactElement} - The component showing the 404 error page.
 */
const NotFound = () => {
    return (
        <Box
            display='flex'
            flexDirection='column'
            justifyContent='center'
            alignItems='center'
            minHeight='80vh' // vertically centers content without taking full screen
            padding={4} // arger padding for better spacing
        >
            <Box
                display='flex'
                flexDirection='column'
                justifyContent='center'
                alignItems='center'
                mb={2}
            >
                <Typography variant='h4' component='div' fontWeight='bold' gutterBottom>
                    The page you are looking for does not exist.
                </Typography>
                {/* <Box
                    component='img'
                    src='/paras_error.png'
                    alt='PARAS Error'
                    sx={{ width: 300 }}
                /> */}
            </Box>
            <Button component={Link} to='/' variant='contained' color='primary'>
                <Typography sx={{ color: 'white.main' }}>
                    Go back to the home page
                </Typography>
            </Button>
        </Box>
    );
};

export default NotFound;