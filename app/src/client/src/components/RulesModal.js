import React from 'react';
import {Box, Typography, Input, Modal, IconButton, Tabs, Tab, Card, CardContent, Stack } from '@mui/material';
import { MdClose } from 'react-icons/md';
import { toast } from 'react-toastify';
import ReactionDrawerContainer from './ReactionDrawer';


const RuleCard = ({ identifier, rule }) => {
        return (
            <Card variant="outlined" sx={{ mb: 2 }}>
                <CardContent>
                    <Typography variant="h6" gutterBottom>
                        {rule.name}
                    </Typography>
                    <Stack direction="column" spacing={2}>
                        <Stack direction="column">
                            <Typography variant="body2" color="text.secondary">
                                Reaction SMARTS:
                            </Typography>
                            <Input
                                value={rule.reaction_smarts}
                                fullWidth
                                multiline
                                readOnly
                            />
                        </Stack>
                        {/* <ReactionDrawerContainer
                            identifier={identifier}
                            // smilesStr={'C=CCBr.[Na+].[I-]>CC(=O)C>C=CCI.[Na+].[Br-]'}
                            smilesStr={rule.reaction_smarts}
                            width={800}
                            height={100}
                        /> */}
                    </Stack>
                </CardContent>
            </Card>
        )
};


const RulesModal = ({
    openRulesModal,
    handleCloseRulesModal,
}) => {
    const [preprocessingRules, setPreprocessingRules] = React.useState([]);
    const [sequencingRules, setSequencingRules] = React.useState([]);
    const [tabValue, setTabValue] = React.useState(0);

    // handle tab change
    const handleTabChange = (event, newValue) => { 
        setTabValue(newValue);
    };

    // handle for compound submission
    const handleGetReactionRules = async () => {
        const data = { };

        try {
            const response = await fetch('/api/get_reaction_rules', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({ data })
            });

            if (!response.ok) {
                throw new Error('Network response was not ok.');
            };

            const json = await response.json();

            if (json.status === 'success') {
                setPreprocessingRules(json.payload.preprocessing_rules.sort((a, b) => a.name.localeCompare(b.name)));
                setSequencingRules(json.payload.sequencing_rules.sort((a, b) => a.name.localeCompare(b.name)));
            } else if (json.status === 'warning') {
                toast.warn(json.message);
            } else if (json.status === 'failure') {
                toast.error(json.message);
            };
        } catch (error) {
            console.error('Error:', error);
            toast.error(error.message);
        };
    };

    // retrieve rules when modal is opened
    React.useEffect(() => {
        if (openRulesModal) {
            handleGetReactionRules();
        };
    }, [openRulesModal]);

    return (
        <Modal open={openRulesModal} onClose={handleCloseRulesModal}>
            <Box
                width='90%'
                bgcolor='white.main'
                mx='auto'
                my={10}
                borderRadius={4}
                boxShadow={3}
                maxHeight='90%'
            >
                <Box
                    sx={{
                        backgroundColor: 'secondary.main',
                        borderTopLeftRadius: '14px',
                        borderTopRightRadius: '14px',
                        display: 'flex',
                        justifyContent: 'space-between',
                        alignItems: 'center',
                        padding: 1,
                    }}
                >
                    <Typography 
                        variant='h5' 
                        gutterBottom
                        sx={{ 
                            color: 'black.main', 
                            textAlign: 'center',
                            pl: 2,
                            pt: 2, 
                        }}
                    >
                        Reaction rule set
                    </Typography>
                    <IconButton onClick={handleCloseRulesModal}>
                        <MdClose size={24} />
                    </IconButton>
                </Box>

                <Box sx={{ width: '100%' }}>
                    <Tabs
                        value={tabValue}
                        onChange={handleTabChange}
                        indicatorColor="primary"
                        textColor="primary"
                        variant="fullWidth"
                        sx={{ borderBottom: 1, borderColor: 'divider' }}
                    >
                        <Tab label="Preprocessing rules" />
                        <Tab label="Sequencing rules" />
                    </Tabs>

                    <Box sx={{ p: 3, overflowY: 'auto', maxHeight: '75vh', mr: 2 }}>
                        {tabValue === 0 && (
                            <Stack direction ="column" spacing={2}>
                                {preprocessingRules.length > 0 ? (
                                    preprocessingRules.map((rule, index) => (
                                        <RuleCard key={index} identifier={`preprocessing-rule-${index}`} rule={rule} />
                                    ))
                                ) : (
                                    <Typography variant="body2" color="text.secondary">
                                        No preprocessing rules available.
                                    </Typography>
                                )}
                            </Stack>
                        )}

                        {tabValue === 1 && (
                            <Stack direction ="column" spacing={2}>
                                {sequencingRules.length > 0 ? (
                                    sequencingRules.map((rule, index) => (
                                        <RuleCard key={index} identifier={`sequencing-rule-${index}`} rule={rule} />
                                    ))
                                ) : (
                                    <Typography variant="body2" color="text.secondary">
                                        No sequencing rules available.
                                    </Typography>
                                )}
                            </Stack>
                        )}
                    </Box>
                </Box>

            </Box>
        </Modal>
    );
};

export default RulesModal;