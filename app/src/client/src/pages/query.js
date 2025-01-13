import React from 'react';
import { Box, Typography } from '@mui/material';
import { useParams } from 'react-router-dom';


const Query = () => {
    // get query from url if any
    const { query } = useParams();

    return (
        <Box
            display='flex'
            flexDirection='column'
            justifyContent='center'
            alignItems='center'

            sx={{
                mt: 5
            }}
        >
            <Typography
                variant='h5'
                align='center'
                gutterBottom
            >
                Querying not yet implemented
            </Typography>

            <Typography
                variant='h6'
                align='center'
                gutterBottom
            >
                query: {query ? query : 'none'}
            </Typography>
        </Box>
    );
};

export default Query;