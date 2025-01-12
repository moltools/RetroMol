import React, { useState } from 'react';
import { Box, Typography, Button, TextField, Grid, Stack, IconButton, Tooltip, Paper, List, ListItem, ListItemButton, ListItemText } from '@mui/material';
import PropTypes from 'prop-types';
import SmilesDrawerContainer from '../components/SmilesDrawer';
import { toast } from 'react-toastify';
import ReplayIcon from '@mui/icons-material/Replay';

const SuccessMessage = ({ message, onReopen, showOnReopen }) => (
    <Box
        display="flex"
        justifyContent="space-between"
        alignItems="center"
        paddingTop={2}
        paddingBottom={2}
        paddingLeft={4}
        paddingRight={4}
        borderRadius={4}
        backgroundColor="#ceccca"
    >
        <Typography variant="body1">
            {message}
        </Typography>
        {showOnReopen && (
            <Tooltip title="Reopen Input">
                <IconButton color="primary" onClick={onReopen} aria-label="Reopen Input">
                    <ReplayIcon />
                </IconButton>
            </Tooltip>
        )}
    </Box>
);

SuccessMessage.propTypes = {
    message: PropTypes.string.isRequired,
    onReopen: PropTypes.func.isRequired,
};

const SmilesInputComponent = ({
    smiles,
    setSmiles,
    handleClear,
    handleSubmitCompound,
    busy,
}) => {
    return (
        <Box padding={4} backgroundColor="#ceccca" borderRadius={4}>
            <Grid container spacing={4} alignItems="flex-start">
                <Grid item xs={12} md={6}>
                    <Stack
                        direction="column"
                        spacing={2}
                    >
                        <TextField
                            label="SMILES"
                            variant="outlined"
                            value={smiles}
                            onChange={(e) => setSmiles(e.target.value)}
                            placeholder="Enter SMILES string"
                            fullWidth
                            aria-label="SMILES input field"
                        />
                        <Stack
                            direction="row"
                            spacing={2}
                        >
                            <Button
                                variant="contained"
                                color="primary"
                                onClick={handleSubmitCompound}
                                disabled={busy}
                            >
                                {busy ? 'Parsing...' : 'Parse'}
                            </Button>
                            <Button
                                variant="contained"
                                color="secondary"
                                onClick={handleClear}
                            >
                                Clear
                            </Button>
                            <Button
                                variant="contained"
                                color="secondary"
                                onClick={() => setSmiles('CCC[C@@H]1C[C@@H](C[C@@H](C[C@@H]2C[C@H](C[C@@H](O2)CC(=O)O1)OC(=O)C=CCCC3=COC(=N3)C=CCNC(=O)OC)C)OC')}
                            >
                                Example: neopeltolide
                            </Button>
                        </Stack>
                    </Stack>
                </Grid>
                <Grid item xs={12} md={6}>
                    <Box
                        display="flex"
                        justifyContent="center"
                        alignItems="center"
                        height="100%"
                        backgroundColor="#fff"
                        borderRadius={1}
                    >
                        <SmilesDrawerContainer
                            identifier="input-compound"
                            smilesStr={smiles}
                            width={300}
                            height={300}
                        />
                    </Box>
                </Grid>
            </Grid>
        </Box>
    );
};

SmilesInputComponent.propTypes = {
    smiles: PropTypes.string.isRequired,
    setSmiles: PropTypes.func.isRequired,
    handleSubmitCompound: PropTypes.func.isRequired,
    busy: PropTypes.bool.isRequired,
};

const ResultComponent = ({ result }) => {
    const smiles = result.smiles;
    const encoding_to_atom_list = result.encoding_to_atom_list;

    // Extract all encoding keys
    const encodingKeys = Object.keys(encoding_to_atom_list);

    // State to manage selected encoding
    const [selectedEncoding, setSelectedEncoding] = useState(
        encodingKeys.length > 0 ? encodingKeys[0] : null
    );

    // Get the corresponding atom list for the selected encoding
    const highlight_atoms = selectedEncoding
        ? encoding_to_atom_list[selectedEncoding].map((atom) => [atom, 'blue'])
        : [];

    // Handler for encoding selection
    const handleSelect = (encoding) => {
        setSelectedEncoding(encoding);
    };

    return (
        <Box padding={4} backgroundColor="#ceccca" borderRadius={4}>
            <Grid container spacing={4} alignItems="flex-start">
                <Grid item xs={12} md={6}>
                    <Paper elevation={3}>
                        <Box padding={2}>
                            <List component="nav" aria-label="encoding options">
                                {encodingKeys.map((encoding) => (
                                    <ListItem key={encoding} disablePadding>
                                        <ListItemButton
                                            selected={selectedEncoding === encoding}
                                            onClick={() => handleSelect(encoding)}
                                        >
                                            {/* <ListItemText primary={encoding} /> */}
                                            <ListItemText primary={result.encoding_to_identity[encoding]} />
                                        </ListItemButton>
                                    </ListItem>
                                ))}
                            </List>
                        </Box>
                    </Paper>
                </Grid>
                <Grid item xs={12} md={6}>
                    <Box
                        display="flex"
                        justifyContent="center"
                        alignItems="center"
                        height="100%"
                        backgroundColor="#fff"
                        borderRadius={1}
                    >
                        <SmilesDrawerContainer
                            identifier="result-compound"
                            smilesStr={smiles}
                            width={300}
                            height={300}
                            highlightAtoms={highlight_atoms}
                        />
                    </Box>
                </Grid>
            </Grid>
        </Box>
    );
}

const Compound = () => {
    const [busy, setBusy] = useState(false);
    const [smiles, setSmiles] = useState('');
    const [result, setResult] = useState(null);

    // collapsed states input components
    const [inputCollapsed, setInputCollapsed] = useState(false);
    const [resultCollapsed, setResultCollapsed] = useState(true);

    // handle for clearing the input field
    const handleClear = () => {
        setSmiles('');
        setInputCollapsed(false);
        setResultCollapsed(true);
        setResult(null);
    };

    // handle for compound submission
    const handleSubmitCompound = async () => {
        setBusy(true);
        setResultCollapsed(true);
        setResult(null);

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
                setInputCollapsed(true);
                setResultCollapsed(false);
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
            <Stack
                direction="column"
                spacing={2}
            >
                <Box>
                    {/* Input container for SMILES string. */}
                    {!inputCollapsed ? (
                        <SmilesInputComponent
                            smiles={smiles}
                            setSmiles={setSmiles}
                            handleClear={handleClear}
                            handleSubmitCompound={handleSubmitCompound}
                            busy={busy} 
                        />
                    ) : (
                        <SuccessMessage
                            message="Compound successfully parsed."
                            onReopen={() => setInputCollapsed(false)}
                            showOnReopen={true}
                        />
                    )}
                </Box>
                
                <Box>
                    {/* Display the result of the compound parsing. */}
                    {!resultCollapsed ? (
                        <ResultComponent result={result} />
                    ) : (
                        <SuccessMessage
                            message={result ? 'Parsing result available for display.' : 'No result to display.'}
                            onReopen={() => setResultCollapsed(false)}
                            showOnReopen={result ? true : false}
                        />
                    )}
                </Box>
            </Stack>

        </Box>
    );
};

export default Compound;