import React from 'react';
import { Box, Typography, Button, Link } from '@mui/material';

/**
 * Home component that displays the home page content.
 * 
 * @returns {React.ReactElement} - The component showing the home page content.
 */
const Home = () => {
    return (
        <Box
            display='flex'
            flexDirection='column'
            justifyContent='center'
            alignItems='center'
            minHeight='80vh' // vertically centers content without taking full screen
            padding={4} // larger padding for better spacing
        >
            <Box display='flex' flexDirection='row' justifyContent='center' alignItems='center'>
                <Box
                    component='img'
                    src='/logo.png'
                    alt='RetroMol logo'
                    sx={{ maxWidth: 150, mb: 2 }}
                />
            </Box>

            <Typography variant='h3' component='div' fontWeight='bold' gutterBottom>
                Welcome to RetroMol!
            </Typography>

            <Typography variant='h6' component='div' color='textSecondary' align='center' gutterBottom>
                Perform cross-modal retrieval between natural product compounds and biosynthetic gene clusters.
            </Typography>
            
            <Box
                sx={{
                    display: 'flex',
                    flexDirection: 'column',
                    justifyContent: 'center',
                    alignItems: 'center',
                    gap: 1.5,
                    mt: 2,
                    mb: 4
                }}
            >
                <Box
                    sx={{
                        textAlign: 'center',
                        display: 'flex',
                        flexDirection: 'row',
                        gap: 2
                    }}
                >
                    <Button
                        variant='contained'
                        color='primary'
                        size='large'
                        href='/session'
                        sx={{ width: 200 }}
                        disabled
                    >
                        <Typography sx={{ color: 'white.main' }}>
                            Retrieve session
                        </Typography>
                    </Button>
                    <Button
                        variant='contained'
                        color='primary'
                        size='large'
                        href='/query'
                        sx={{ width: 200 }}
                    >
                        <Typography sx={{ color: 'white.main'}}>
                            manual query
                        </Typography>
                    </Button>
                </Box>
                <Box
                    sx={{
                        textAlign: 'center',
                        display: 'flex',
                        flexDirection: 'row',
                        gap: 2
                    }}
                >
                    <Button
                        variant='contained'
                        color='primary'
                        size='large'
                        href='/compound'
                        sx={{ width: 200 }}
                    >
                        <Typography sx={{ color: 'white.main' }}>
                            Parse compound
                        </Typography>
                    </Button>
                    <Button
                        variant='contained'
                        color='primary'
                        size='large'
                        href='/antismash'
                        sx={{ width: 200 }}
                        disabled
                    >
                        <Typography sx={{ color: 'white.main'}}>
                            Parse antiSMASH
                        </Typography>
                    </Button>
                </Box>
            </Box>

            <Typography variant='body1' align='center' color='textSecondary' gutterBottom>
                Want to learn more about the research behind RetroMol?
                <Link 
                    href='/publication' 
                    underline='hover' 
                    sx={{ marginLeft: 1, fontWeight: 'bold' }}
                >
                    Read our publication.
                </Link>
            </Typography>

            <Typography variant='body1' align='center' color='textSecondary'>
                Have something to contribute?
                <Link 
                    href='https://github.com/moltools/RetroMol'
                    target='_blank'
                    underline='hover'
                    sx={{ marginLeft: 1, fontWeight: 'bold' }}
                >
                    Visit our GitHub page.
                </Link>
            </Typography>
        </Box>
    );
};

export default Home;