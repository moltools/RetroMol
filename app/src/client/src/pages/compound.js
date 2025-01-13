import React, { useState } from 'react';
import { Box, Typography, Button, TextField, Grid, Stack, IconButton, Tooltip, Paper, List, ListItem, ListItemButton, ListItemText } from '@mui/material';
import PropTypes from 'prop-types';
import SmilesDrawerContainer from '../components/SmilesDrawer';
import { toast } from 'react-toastify';
import ReplayIcon from '@mui/icons-material/Replay';
import ExitIcon from '@mui/icons-material/ExitToApp';

const SuccessMessage = ({ message, onReopen, showOnReopen }) => (
    <Box
        sx={{
            backgroundColor:"#ceccca",
            borderRadius: 4,
            paddingTop: 2,
            paddingBottom: 2,
            paddingLeft: 4,
            paddingRight: 4,
            alignItems: 'center',
            justifyContent: 'space-between',
            display: 'flex',
        }}
    >
        <Typography variant="body1">
            {message}
        </Typography>
        {showOnReopen && (
            <Tooltip title="Reopen input form." arrow>
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
        <Box padding={4} sx={{backgroundColor: "#ceccca", borderRadius: 4}}>
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
                                <Typography variant="body1" sx={{ color: 'common.white' }}>
                                    {busy ? 'Parsing...' : 'Parse'}
                                </Typography>
                            </Button>
                            <Button
                                variant="contained"
                                color="secondary"
                                onClick={handleClear}
                            >
                                <Typography variant="body1" sx={{ color: 'common.black' }}>
                                    Clear
                                </Typography>
                            </Button>
                            <Button
                                variant="contained"
                                color="secondary"
                                onClick={() => setSmiles("CCC[C@@H]1C[C@@H](C[C@@H](C[C@@H]2C[C@H](C[C@@H](O2)CC(=O)O1)OC(=O)/C=C\\CCC3=COC(=N3)/C=C\\CNC(=O)OC)C)OC")}
                            >
                                <Typography variant="body1" sx={{ color: 'common.black' }}>
                                    Example: neopeltolide
                                </Typography>
                            </Button>
                        </Stack>
                    </Stack>
                </Grid>
                <Grid item xs={12} md={6}>
                    <Paper
                        elevation={3}
                        sx={{
                            backgroundColor: "#fff",
                            display: 'flex',
                            justifyContent: 'center',
                            alignItems: 'center',
                            borderRadius: 1
                        }}
                    >
                        <SmilesDrawerContainer
                            identifier="input-compound"
                            smilesStr={smiles}
                            width={500}
                            height={500}
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
        // ? encoding_to_atom_list[selectedEncoding].map((atom) => [atom, '#3d7dca'])
        ? encoding_to_atom_list[selectedEncoding].map((atom) => [atom, '#ceccca'])
        : [];

    // const linear_highlight_atoms = selectedEncoding
    //     ? result.encoding_to_highlights[selectedEncoding]
    //     : [];

    const [linearHighlightsAtoms, setLinearHighlightsAtoms] = useState([]);
    const [selectedMotifIndex, setSelectedMotifIndex] = useState(null);

    // Handler for encoding selection
    const handleSelect = (encoding) => {
        setSelectedEncoding(encoding);
        setForwardedEncoding(encoding);
        setLinearSmiles(result.encoding_to_smiles[encoding]);
        setSelectedMotifIndex(null);
        setLinearHighlightsAtoms([]);
    };

    return (
        <Box padding={4} sx={{backgroundColor: "#ceccca", borderRadius: 4}}> 
            <Grid container spacing={4} alignItems="flex-start">
                <Grid item xs={12} md={6}>
                    {/* <Typography variant="body1" sx={{ height: 50, marginLeft: 1}}>
                        Select backbone:
                    </Typography> */}
                    <Stack
                        direction="column"
                        spacing={2}
                    >
                        <Paper elevation={3} sx ={{ padding: 0 }}>
                            <Box padding={0} sx={{ overflow: 'auto' }}>
                                <Box
                                    sx={{
                                        borderTopLeftRadius: 4,
                                        borderTopRightRadius: 4,
                                        backgroundColor: 'primary.main',
                                        color: 'common.white',
                                    }}
                                >
                                    <Box padding={2} sx={{ fontWeight: 'regular', color: 'inherit', textAlign: 'center' }}>
                                        SELECT PRIMARY SEQUENCE
                                    </Box>
                                </Box>
                                <List component="nav" aria-label="encoding options">
                                    {encodingKeys.map((encoding, index) => (
                                        <ListItem key={encoding} disablePadding>
                                            <ListItemButton
                                                selected={selectedEncoding === encoding}
                                                onClick={() => handleSelect(encoding)}
                                            >
                                                {/* <ListItemText primary={result.encoding_to_identity[encoding]} /> */}
                                                <ListItemText primary={`Primary sequence ${index + 1}`} />
                                            </ListItemButton>
                                        </ListItem>
                                    ))}
                                </List>
                            </Box>
                        </Paper>
                        <Tooltip
                            title="Forward primary sequence to query page."
                            arrow
                        >
                            <Button
                                variant="contained"
                                color="secondary"
                                sx={{
                                    width: '100%',
                                    height: 50,
                                    center: 'center',
                                    alignItems: 'center',
                                }}
                                onClick={() => {
                                    // open /query/:query? page with selected primary sequence
                                    const query = result.encoding_to_primary_sequence[selectedEncoding].map((motif) => motif.name).join('>');
                                    window.open(`/query/${query}`, '_blank');
                                }}
                            >
                                <Typography variant="body1" sx={{ color: 'common.black', fontWeight: 'bold' }}>
                                        Forward primary sequence
                                </Typography>
                                <ExitIcon sx={{ marginLeft: 1 }} />
                            </Button>
                        </Tooltip>
                    </Stack>
                </Grid>
                <Grid item xs={12} md={6}>
                    <Paper
                        sx={{backgroundColor:"#fff", borderRadius: 1}}  
                        elevation={3}
                    >
                        <Stack direction="column">
                            {/* box on top has two tabs left and right */}
                            <Box sx={{ height: 50, width: '100%' }}>
                                <Stack
                                    direction="row"
                                    spacing={0}
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
                                            boxShadow: molviewTab === 'overview' ? 'inset 0px 2px 5px rgba(0,0,0,1)' : undefined,
                                            transition: 'all 0.3s ease-in-out',
                                            '&:active': {
                                                boxShadow: molviewTab === 'overview' ? 'inset 0px 2px 5px rgba(0,0,0,1)' : undefined,
                                                transform: molviewTab === 'overview' ? 'translateY(2px)' : 'none',
                                            },
                                        }}
                                    >
                                        <Typography variant="body1" sx={{ color: 'common.white' }}>
                                            Overview
                                        </Typography>
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
                                            boxShadow: molviewTab === 'linear' ? 'inset 0px 2px 5px rgba(0,0,0,1)' : undefined,
                                            transition: 'all 0.2s ease-in-out',
                                            '&:active': {
                                                boxShadow: molviewTab === 'linear' ? 'inset 0px 2px 5px rgba(0,0,0,1)' : undefined,
                                                transform: molviewTab === 'linear' ? 'translateY(2px)' : 'none',
                                            },
                                        }}
                                    >
                                        <Typography variant="body1" sx={{ color: 'common.white' }}>
                                            Linear view
                                        </Typography>
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
                                        width={500}
                                        height={500}
                                        highlightAtoms={highlight_atoms}
                                    />
                                ) : (
                                    <SmilesDrawerContainer
                                        identifier="result-compound"
                                        smilesStr={linearSmiles}
                                        width={500}
                                        height={500}
                                        highlightAtoms={linearHighlightsAtoms}
                                    />
                                )}
                            </Box>
                            <Box sx={{ height: 90, paddingBottom: 1.5 }}>
                                <Stack direction="column" spacing={0.1}>
                                    <Typography variant="body1" sx={{ height: 25, paddingLeft: 2.5, paddingBottom: 1 }}>
                                        Primary sequence:
                                    </Typography>
                                    <Stack
                                        direction="row" 
                                        spacing={0.3}
                                        sx={{
                                            display: 'flex',
                                            justifyContent: 'flex-start',
                                            alignItems: 'flex-start',
                                            height: '100%',
                                            overflowX: 'auto',
                                            overflowY: 'hidden',
                                            flexWrap: 'nowrap',
                                            paddingLeft: 1,
                                            paddingBottom: 2.5,
                                            paddingRight: 1,
                                        }}
                                    >
                                        {result.encoding_to_primary_sequence[selectedEncoding].map((motif, index) => (
                                            <Tooltip key={index} title="Highlight motif in linear view." arrow enterDelay={500}>
                                            <Button 
                                                key={index} 
                                                sx={{ 
                                                    padding: 1,
                                                    backgroundColor: 'secondary.main',
                                                    flexShrink: 0,
                                                    boxShadow: 1,
                                                    borderRadius: 2,
                                                    boxShadow: selectedMotifIndex === index ? 'inset 0px 2px 5px rgba(0,0,0,1)' : undefined,
                                                    transition: 'all 0.2s ease-in-out',
                                                    '&:active': {
                                                        boxShadow: selectedMotifIndex === index ? 'inset 0px 2px 5px rgba(0,0,0,1)' : undefined,
                                                        transform: selectedMotifIndex === index ? 'translateY(2px)' : 'none',
                                                    },
                                                }}
                                                onClick={() => {
                                                    // if it is already selected, deselect it
                                                    if (selectedMotifIndex === index) {
                                                        setLinearHighlightsAtoms([]);
                                                        setSelectedMotifIndex(null);
                                                        return;
                                                    }

                                                    const atomsToHighlight = result.encoding_to_highlights[selectedEncoding][index];
                                                    const newAtomsToHighlight = atomsToHighlight.map((atom) => [atom, '#ceccca']);
                                                    setLinearHighlightsAtoms(newAtomsToHighlight);
                                                    setSelectedMotifIndex(index);
                                                }}
                                            >
                                                <Typography 
                                                    variant="body1" 
                                                    noWrap
                                                    sx={{ fontWeight: 'bold' }}
                                                >
                                                    {motif.name}
                                                </Typography>
                                            </Button>
                                            </Tooltip>
                                        ))}
                                    </Stack>
                                </Stack>
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
                // check if any json.payload.encoding_to_identity value is equal to _backbone
                const backbone = Object.values(json.payload.encoding_to_identity).includes('_backbone');
                if (!backbone) {
                    toast.warn('No backbone found in the compound. If the compound is a modular natural product, the rule set is likely insufficient to parse this compound.');
                    setBusy(false);
                    return;
                }
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
                            showOnReopen={false}
                        />
                    )}
                </Box>
                    
            </Stack>

        </Box>
    );
};

export default Compound;