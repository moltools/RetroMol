import React, { useState } from 'react';
import { Box, Typography, Button, Link } from '@mui/material';
import SmilesDrawerContainer from '../components/SmilesDrawer';
import { toast } from 'react-toastify';

/**
 * Compound component that shows the workspace for parsing compounds.
 * 
 * @returns {React.ReactElement} - The component showing the compound page content.
 */
const Compound = () => {
    const [busy, setBusy] = useState(false);
    const [smiles, setSmiles] = useState('');
    const [result, setResult] = useState({});

    // handle for compound submission
    const handleSubmitCompound = async () => {
        setBusy(true);

        const data = { smiles: smiles };

        try {
            const response = await fetch('/api/submit_compound', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({ data })
            });

            if (!response.ok) {
                throw new Error('Network response was not ok.');
            };

            const json = await response.json();

            if (json.status === 'success') {
                toast.success(json.message);
                setResult(json.payload);
            } else if (json.status === 'warning') {
                toast.warn(json.message);
            } else if (json.status === 'failure') {
                toast.error(json.message);
            };
        } catch (error) {
            console.error('Error:', error);
            toast.error(error.message);
        };

        setBusy(false);
    };

    return (
        <Box
            display='flex'
            flexDirection='column'
            justifyContent='left'
            padding={4}
            margin='auto'
        >
            {/* Input container for SMILES string. */}
            <Box
                display='flex'
                flexDirection='row'
                justifyContent='left'
                alignItems='center'
                padding={2}
            >
                <Typography variant='h6'>
                    SMILES:
                </Typography>
                <input
                    type='text'
                    value={smiles}
                    onChange={(e) => setSmiles(e.target.value)}
                    style={{ width: '100%' }}
                />

                {/* Button to parse SMILES string. */}
                <Button
                    variant='contained'
                    color='primary'
                    style={{ marginLeft: '1rem' }}
                    onClick={handleSubmitCompound}
                    disabled={busy}
                >
                    Parse
                </Button>

                {/* Button to clear the input. */}
                <Button
                    variant='contained'
                    color='secondary'
                    style={{ marginLeft: '1rem' }}
                    onClick={() => setSmiles('')}
                >
                    Clear
                </Button>

                {/* Draw SMILES when typing. */}
                <SmilesDrawerContainer 
                    identifier={'input-compound'}
                    smilesStr={smiles} 
                    width={300} 
                    height={100}
                />
            </Box>

            {/* Display result in text format. */}
            <Box
                display='flex'
                flexDirection='column'
                justifyContent='left'
                alignItems='left'
                padding={2}
            >
                <Typography variant='h6'>
                    Result:
                </Typography>
                <pre>
                    {JSON.stringify(result, null, 2)}
                </pre> 
            </Box>

        </Box>
    );
};

export default Compound;