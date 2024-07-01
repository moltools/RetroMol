import React, { useState, useEffect } from "react";
import { toast } from "react-toastify";
import { Link } from "react-router-dom";
import { 
    AppBar, 
    Box, 
    Button,
    Container,
    Divider,
    Drawer, 
    Grid,
    IconButton, 
    List, 
    ListItem, 
    ListItemIcon, 
    ListItemText, 
    Modal,
    ToggleButton,
    ToggleButtonGroup,
    Toolbar, 
    Typography 
} from "@mui/material";
import { 
    BugReport as BugReportIcon,
    Home as HomeIcon, 
    Info as InfoIcon, 
    Menu as MenuIcon,
    Close as CloseIcon
} from "@mui/icons-material";

import Alignment from "../components/dashboard/Alignment";
import InputForm from "../components/dashboard/InputForm";
import LoadingOverlay from "../components/common/LoadingOverlay";
import Page from "../components/common/Page";
import QueryForm from "../components/dashboard/QueryForm";
import ResultsDisplay from "../components/dashboard/ResultsDisplay";
import Status from "../components/common/Status";

const SidebarButtonHome = () => {
    return (
        <ListItem 
            type="button" 
            component={Link} 
            to="/"
            sx={{ textDecoration: "none", color: "#222" }}
        >
            <ListItemIcon>
                <HomeIcon sx={{ color: "#222" }} />
            </ListItemIcon>
            <ListItemText primary="Home" />
        </ListItem>
    );
};

const SidebarButtonAbout = () => {
    return (
        <ListItem 
            type="button" 
            component={Link} 
            to="/about"
            sx={{ textDecoration: "none", color: "#222" }}
        >
            <ListItemIcon>
                <InfoIcon sx={{ color: "#222" }} />
            </ListItemIcon>
            <ListItemText primary="About" />
        </ListItem>
    );
};

const SidebarButtonBugReport = () => {
    return (
        <ListItem 
            type="button" 
            component="a" 
            href="https://github.com/moltools/RetroMol/issues" 
            target="_blank"
            sx={{ textDecoration: "none", color: "#222" }}
        >
            <ListItemIcon>
                <BugReportIcon sx={{ color: "#222" }} />
            </ListItemIcon>
            <ListItemText primary="Report Bug" />
        </ListItem>
    );
};

const RetroMol = () => {
    const [isBusy, setIsBusy] = useState(false);
    const [queryOnly, setQueryOnly] = useState(false);
    const [isDrawerOpen, setIsDrawerOpen] = useState(false);
    const [version, setVersion] = useState("0.0.0");
    const [statusServer, setStatusServer] = useState(true);
    const [statusDatabase, setStatusDatabase] = useState(false);
    const [parsedResults, setParsedResults] = useState([]);
    const [selectedResultIndex, setSelectedResultIndex] = useState(null);
    
    const [query, setQuery] = useState([]);
    const [resultModalOpen, setResultModalOpen] = useState(false);
    const [queryResult, setQueryResult] = useState({});

    const toggleDrawer = () => {
        setIsDrawerOpen(!isDrawerOpen);
    };

    const fetchVersion = async () => {
        try {
            const response = await fetch("/api/fetch_version");
            const data = await response.json();
            if (data.status === "success") {
                setVersion(data.payload.version);
            };
        } catch (error) {
            console.error(error);
        };
    };

    const checkStatusServer = async () => {
        try {
            const response = await fetch("/api/fetch_server_status");
            const data = await response.json();
            if (data.status === "success") {
                setStatusServer(true);
            } else {
                setStatusServer(false);
            };
        } catch (error) {
            setStatusServer(false);
            console.error(error);
        };
    };

    const checkStatusDatabase = async () => {
        try {
            const response = await fetch("/api/fetch_database_status");
            const data = await response.json();
            if (data.status === "success") {
                setStatusDatabase(true);
            } else {
                setStatusDatabase(false);
            };
        } catch (error) {
            setStatusDatabase(false);
            console.error(error);
        };
    };

    useEffect(() => {
        fetchVersion();

        // toast.warning(
        //     "SMILES parsing is under active development and may not work as expected. Please report any issues.",
        //     {
        //         position: "top-center",
        //         hideProgressBar: false,
        //         closeOnClick: true,
        //         pauseOnHover: true,
        //         draggable: true,
        //         progress: undefined,
        //         autoClose: false,
        //     }
        // );

        if (isDrawerOpen) {
            checkStatusServer();
            checkStatusDatabase();
        }
    }, [isDrawerOpen]);

    const submit = async (data) => {
        setIsBusy(true);

        // reset selectedResultIndex
        setParsedResults([]);
        setSelectedResultIndex(null);

        try {
            const response = await fetch("/api/parse_submission", {
                method: "POST",
                headers: {"Content-Type": "application/json"},
                body: JSON.stringify(data)
            });

            const result = await response.json();
            if (result.status === "success") {
                toast.success(result.message);
                setParsedResults(result.payload.queries);
            } else if (result.status === "warning") {
                toast.warn(result.message);
            } else {
                toast.error(result.message);
            };
            
        } catch (error) {
            console.error(error);
            toast.error("Something went wrong. Please file a bug report.");
        };

        setIsBusy(false);
    };

    const submitQuery = async (data) => {
        setIsBusy(true);

        // setParsedResults([]);

        try {
            const response = await fetch("/api/query_submission", {
                method: "POST",
                headers: {"Content-Type": "application/json"},
                body: JSON.stringify(data)
            });

            const result = await response.json();
            if (result.status === "success") {
                toast.success(result.message);
                setQueryResult(result.payload["motif_codes"]);
                setResultModalOpen(true);
            } else if (result.status === "warning") {
                toast.warn(result.message);
            } else {
                toast.error(result.message);
            };
            
        } catch (error) {
            console.error(error);
            toast.error("Something went wrong. Please file a bug report.");
        };

        setIsBusy(false);
    };

    const handleDownload = () => {
        const filename = "results.json";
        const dataJson = JSON.stringify(queryResult, null, 4);
        const blob = new Blob([dataJson], { type: "application/json" });
        const url = URL.createObjectURL(blob);
        const link = document.createElement("a");
        link.href = url;
        link.download = filename;
        link.click();
    };

    return (
        <Page sx={{ display: "flex", width: "100%", flexDirection: "column" }}>
            <AppBar position="fixed">
                <Toolbar sx={{ display: "flex", justifyContent: "space-between" }}>
                    <Box sx={{ display: "flex", justifyContent: "left", width: "100%", alignItems: "center" }}>
                        <IconButton
                            color="inherit"
                            aria-label="open drawer"
                            edge="start"
                            onClick={toggleDrawer}
                            sx={{ mr: 2 }}
                        >
                            <MenuIcon />
                        </IconButton>
                        <Box sx={{ display: "flex", alignItems: "center" }}>
                            <Box sx={{ display: "flex", alignItems: "center", mr: 2 }}>
                                <img src="/retromol.svg" alt="RetroMol" style={{ width: 40, height: 40 }} />
                            </Box>
                            <Typography variant="h6" noWrap>
                                {`RetroMol (${version})`}
                            </Typography>
                        </Box>
                    </Box>
                </Toolbar>
            </AppBar>
            <Drawer
                variant="temporary"
                open={isDrawerOpen}
                onClose={toggleDrawer}
                ModalProps={{ keepMounted: true }}
                sx={{ '& .MuiDrawer-paper': { width: 240 } }}
            >
                <List>
                    <SidebarButtonHome />
                    <SidebarButtonAbout />
                    <SidebarButtonBugReport />
                    <Divider sx={{ mt: 2, mb: 2 }} />
                    <Status statusName="Server" status={statusServer} />
                    <Status statusName="Database" status={statusDatabase} />
                </List>
            </Drawer>
            <Box
                component="main"
                sx={{ flexGrow: 1, p: 3 }}
            >
                <Toolbar />
                <Modal
                    open={resultModalOpen}
                    onClose={() => setResultModalOpen(false)}
                    aria-labelledby="modal-modal-title"
                    aria-describedby="modal-modal-description"
                >
                    <Box
                        sx={{
                            position: "absolute",
                            top: "50%",
                            left: "50%",
                            transform: "translate(-50%, -50%)",
                            width: "90%",
                            height: "80%",
                            bgcolor: "background.paper",
                            boxShadow: 24,
                            borderRadius: 4,
                        }}
                    >   
                        <Box sx={{ 
                            display: "flex", 
                            justifyContent: "space-between", 
                            backgroundColor: "#555", 
                            color: "white", 
                            p: 2,
                            borderTopLeftRadius: 16, 
                            borderTopRightRadius: 16
                        }}>
                            <Typography id="modal-modal-title" variant="h6" component="h2">
                                Query Result
                            </Typography>
                            <Button
                                variant="contained"
                                onClick={handleDownload}
                            >
                                Download results
                            </Button>
                        </Box>
                        <Box 
                            sx={{ 
                                p: 2, 
                                overflow: "auto", 
                                height: "calc(100% - 64px)" 
                            }}
                        >
                            {queryResult.length === 0 ? (
                                <Typography id="modal-modal-description" sx={{ mt: 2 }}>
                                    No results found.
                                </Typography>
                            ) : (
                                <Alignment data={queryResult} />
                            )}
                        </Box>
                    </Box>
                </Modal>
                <Container>
                    {isBusy && <LoadingOverlay />}
                    <Grid container spacing={3}>
                        <Grid item xs={12}>
                            <ToggleButtonGroup
                                value={queryOnly}
                                exclusive
                                onChange={() => setQueryOnly(!queryOnly)}
                                aria-label="query only"
                                sx={{
                                    borderRadius: 2,
                                    boxShadow: 3,
                                }}
                            >
                                <ToggleButton value={false} aria-label="query and parse">
                                    Parse and query
                                </ToggleButton>
                                <ToggleButton value={true} aria-label="query only">
                                    Query only
                                </ToggleButton>
                            </ToggleButtonGroup>
                        </Grid>
                        {!queryOnly && (
                            <Grid item xs={12}>
                                <InputForm onSubmit={(data) => submit(data)} />
                            </Grid>
                        )}
                        {!queryOnly && (
                            <Grid item xs={12}>
                                <ResultsDisplay 
                                    results={parsedResults} 
                                    setResults={setParsedResults} 
                                    selectedResultIndex={selectedResultIndex}
                                    setSelectedResultIndex={setSelectedResultIndex}
                                />
                            </Grid>
                        )}
                        <Grid item xs={12}>
                            <QueryForm
                                results={parsedResults}
                                selectedResultIndex={selectedResultIndex}
                                setSelectedResultIndex={setSelectedResultIndex}
                                columns={query}
                                setColumns={setQuery}
                                submit={submitQuery}
                            />
                        </Grid>
                    </Grid>
                </Container>
            </Box>
        </Page>
    );
};

export default RetroMol;