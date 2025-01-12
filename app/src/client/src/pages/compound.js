import React, { useState } from 'react';
import { Box, Typography, Button, TextField, Grid, Stack, IconButton, Tooltip } from '@mui/material';
import PropTypes from 'prop-types';
import SmilesDrawerContainer from '../components/SmilesDrawer';
import { toast } from 'react-toastify';
import ReplayIcon from '@mui/icons-material/Replay';
import { Graph, VF2Matcher } from '../components/VF2Matcher';

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
    return (
        <Box padding={4} backgroundColor="#ceccca" borderRadius={4}>
            <Grid container spacing={4} alignItems="flex-start">
                <Grid item xs={12} md={6}>
                    <Stack
                        direction="column"
                        spacing={2}
                    >
                        <Typography variant="body1">
                            {JSON.stringify(result.encoding_to_identity, null, 2)}
                        </Typography>
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
                            smilesStr={result.smiles}
                            width={300}
                            height={300}
                        />
                    </Box>
                </Grid>
            </Grid>
        </Box>
    );
}

function createMolecule(bonds) {
    const graph = new Graph();
    // Add all bonds to the graph
    bonds.forEach(bond => graph.addEdge(bond.source, bond.target, bond.bond));
    return graph;
}

const Compound = () => {

    // Define atoms and bonds for the first ethanol molecule
    const bonds1 = [
        { source: 'C1', target: 'C2', bond: 'single' },
        { source: 'C2', target: 'O1', bond: 'single' },
        { source: 'C1', target: 'H1', bond: 'single' },
        { source: 'C1', target: 'H2', bond: 'single' },
        { source: 'C1', target: 'H3', bond: 'single' },
        { source: 'C2', target: 'H4', bond: 'single' },
        { source: 'C2', target: 'H5', bond: 'single' },
        { source: 'O1', target: 'H6', bond: 'single' }
    ];

    // Define atoms and bonds for the second ethanol molecule (different indexing)
    const bonds2 = [
        { source: 'A', target: 'B', bond: 'single' },
        { source: 'B', target: 'C', bond: 'single' },
        { source: 'A', target: 'D', bond: 'single' },
        { source: 'A', target: 'E', bond: 'single' },
        { source: 'A', target: 'F', bond: 'single' },
        { source: 'B', target: 'G', bond: 'single' },
        { source: 'B', target: 'H', bond: 'single' },
        { source: 'C', target: 'I', bond: 'single' }
    ];

    // Create molecular graphs
    const mol1 = createMolecule(bonds1);
    const mol2 = createMolecule(bonds2);

    // Initialize VF2 Matcher
    const matcher = new VF2Matcher(mol1, mol2);

    // Check isomorphism
    if (matcher.isIsomorphic()) {
        console.log("Graphs are isomorphic.");
        console.log("Mapping:", matcher.getMapping());
    } else {
        console.log("Graphs are not isomorphic.");
    }

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