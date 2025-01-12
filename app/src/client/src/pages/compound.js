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
            <Tooltip title="Reopen input box" arrow>
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
                    <Paper
                        backgroundColor="#fff"
                        borderRadius={1}
                        elevation={3}
                        sx={{
                            display: 'flex',
                            justifyContent: 'center',
                            alignItems: 'center',
                        }}
                    >
                        <SmilesDrawerContainer
                            identifier="input-compound"
                            smilesStr={smiles}
                            width={300}
                            height={300}
                        />
                    </Paper>
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

const ResultComponent = ({ setForwardedEncoding, result }) => {
    const smiles = result.smiles;
    const encoding_to_atom_list = result.encoding_to_atom_list;

    const [molviewTab, setMolviewTab] = useState('overview');

    // Extract all encoding keys, but only if 'result.encoding_to_identity[encoding]' is equal to _backbone
    const encodingKeys= Object.keys(encoding_to_atom_list).filter((encoding) => result.encoding_to_identity[encoding] === '_backbone');

    // State to manage selected encoding
    const [selectedEncoding, setSelectedEncoding] = useState(
        encodingKeys.length > 0 ? encodingKeys[0] : null
    );

    const [linearSmiles, setLinearSmiles] = useState(
        selectedEncoding ? result.encoding_to_smiles[selectedEncoding] : null
    );

    // Get the corresponding atom list for the selected encoding
    const highlight_atoms = selectedEncoding
        ? encoding_to_atom_list[selectedEncoding].map((atom) => [atom, '#3d7dca'])
        : [];

    const linear_highlight_atoms = selectedEncoding
        ? result.encoding_to_highlights[selectedEncoding]
        : [];

    // Handler for encoding selection
    const handleSelect = (encoding) => {
        setSelectedEncoding(encoding);
        setForwardedEncoding(encoding);
        setLinearSmiles(result.encoding_to_smiles[encoding]);
    };

    return (
        <Box padding={4} backgroundColor="#ceccca" borderRadius={4}>
            <Grid container spacing={4} alignItems="flex-start">
                <Grid item xs={12} md={6}>
                    <Typography variant="body1" sx={{ height: 50, marginLeft: 1}}>
                        Select backbone:
                    </Typography>
                    <Paper 
                        elevation={3}
                        sx={{
                            height: 300,
                            overflow: 'auto',
                        }}
                    >
                        <Box padding={2}>
                            <List component="nav" aria-label="encoding options">
                                {encodingKeys.map((encoding, index) => (
                                    <ListItem key={encoding} disablePadding>
                                        <ListItemButton
                                            selected={selectedEncoding === encoding}
                                            onClick={() => handleSelect(encoding)}
                                        >
                                            {/* <ListItemText primary={result.encoding_to_identity[encoding]} /> */}
                                            <ListItemText primary={`Backbone ${index + 1}`} />
                                        </ListItemButton>
                                    </ListItem>
                                ))}
                            </List>
                        </Box>
                    </Paper>
                </Grid>
                <Grid item xs={12} md={6}>
                    <Paper
                        backgroundColor="#fff"
                        borderRadius={1}
                        elevation={3}
                    >
                        <Stack direction="column">
                            {/* box on top has two tabs left and right */}
                            <Box sx={{ height: 50, width: '100%' }}>
                                <Stack
                                    direction="row"
                                    spacing={0.1}
                                    sx={{
                                        display: 'flex',
                                        justifyContent: 'center',
                                        alignItems: 'center',
                                        height: '100%',
                                        width: '100%',
                                    }}
                                >
                                    <Button
                                        variant="contained"
                                        color="primary"
                                        onClick={() => setMolviewTab('overview')}
                                        sx={{
                                            borderTopLeftRadius: '4px',
                                            borderTopRightRadius: 0,
                                            borderBottomLeftRadius: 0,
                                            borderBottomRightRadius: 0,
                                            flex: 1,
                                            height: '50px',
                                        }}
                                    >
                                        Overview
                                    </Button>
                                    <Button
                                        variant="contained"
                                        color="primary"
                                        onClick={() => setMolviewTab('linear')}
                                        sx={{
                                            borderTopLeftRadius: 0,
                                            borderTopRightRadius: '4px',
                                            borderBottomLeftRadius: 0,
                                            borderBottomRightRadius: 0,
                                            flex: 1,
                                            height: '50px',
                                        }}
                                    >
                                        Backbone 
                                    </Button>
                                </Stack>
                            </Box>
                            <Box
                                sx={{
                                    display: 'flex',
                                    justifyContent: 'center',
                                    alignItems: 'center',
                                }}
                            >
                                {molviewTab === 'overview' ? (
                                    <SmilesDrawerContainer
                                        identifier="result-compound"
                                        smilesStr={smiles}
                                        width={300}
                                        height={300}
                                        highlightAtoms={highlight_atoms}
                                    />
                                ) : (
                                    <SmilesDrawerContainer
                                        identifier="result-compound"
                                        smilesStr={linearSmiles}
                                        width={300}
                                        height={300}
                                        highlightAtoms={linear_highlight_atoms}
                                    />
                                )}
                            </Box>
                        </Stack>
                    </Paper>  
                </Grid>
            </Grid>
        </Box>
    );
}

ResultComponent.propTypes = {
    setForwardedEncoding: PropTypes.func.isRequired,
    result: PropTypes.object.isRequired,
};

const Compound = () => {
    const [busy, setBusy] = useState(false);
    const [smiles, setSmiles] = useState('');
    const [result, setResult] = useState(null);
    const [encoding, setEncoding] = useState(null);

    // collapsed states input components
    const [inputCollapsed, setInputCollapsed] = useState(false);
    const [resultCollapsed, setResultCollapsed] = useState(true);
    const [encodingCollapsed, setEncodingCollapsed] = useState(true);

    // handle for clearing the input field
    const handleClear = () => {
        setSmiles('');
        setInputCollapsed(false);
        setResultCollapsed(true);
        setResult(null);
        setEncodingCollapsed(true);
        setEncoding(null);
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
                        <ResultComponent 
                            setForwardedEncoding={setEncoding}
                            result={result} 
                        />
                    ) : (
                        <SuccessMessage
                            message={result ? 'Parsing result available for display.' : 'No result to display.'}
                            onReopen={() => setResultCollapsed(false)}
                            showOnReopen={true}
                        />
                    )}
                </Box>
                    
            </Stack>

        </Box>
    );
};

export default Compound;