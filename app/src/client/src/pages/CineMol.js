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
    IconButton, 
    List, 
    ListItem, 
    ListItemIcon, 
    ListItemText, 
    Toolbar, 
    Typography,
    Slider,
    CssBaseline
} from "@mui/material";
import { 
    BugReport as BugReportIcon,
    Home as HomeIcon, 
    Info as InfoIcon, 
    Menu as MenuIcon,
    Close as CloseIcon,
    DragHandle as DragHandleIcon,
    SwitchRight
} from "@mui/icons-material";
import Draggable from 'react-draggable';

import LoadingOverlay from "../components/common/LoadingOverlay";
import Page from "../components/common/Page";
import Status from "../components/common/Status";

// Example SDF string for penicillin G.
const exampleSdfString = `5904
  -OEChem-12182318453D

 41 43  0     1  0  0  0  0  0999 V2000
   -0.8019    1.2308    0.5170 S   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2842   -2.5451   -1.2026 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.3517    1.0760   -0.8170 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.1157   -0.6970    0.5961 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.1598   -2.0405    1.2167 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4781   -0.7369    0.3018 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.5677   -1.3807   -0.3375 N   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3100   -0.4177    1.1279 C   0  0  1  0  0  0  0  0  0  0  0  0
   -2.5670    1.6679    0.1134 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1671    0.3360   -0.3862 C   0  0  2  0  0  0  0  0  0  0  0  0
   -0.6142   -1.6325    0.4990 C   0  0  1  0  0  0  0  0  0  0  0  0
   -1.9181   -1.8261   -0.3007 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.2193    2.2224    1.3871 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5771    2.7313   -0.9863 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.6296    0.1576   -0.1297 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.8730   -1.6107    0.1024 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.9020   -1.2655   -0.9563 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.8679   -0.1949   -0.5225 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.5670    1.1373   -0.7687 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.0457   -0.5556    0.1171 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.4619    2.1290   -0.3670 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.9405    0.4362    0.5190 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.6486    1.7786    0.2769 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5167   -0.4852    2.1999 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.9961    0.1926   -1.4616 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4509   -2.4585    1.2025 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.2655    1.4761    2.1882 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6724    3.0921    1.7709 H   0  0  0  0  0  0  0  0  0  0  0  0
   -4.2446    2.5572    1.1944 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0666    3.6467   -0.6644 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0842    2.3779   -1.8995 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.6008    3.0119   -1.2561 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.4359   -1.0113   -1.2762 H   0  0  0  0  0  0  0  0  0  0  0  0
   -6.3150    0.9762   -0.6611 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.4453   -2.1891   -1.1950 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.4229   -0.9646   -1.8960 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.6443    1.4206   -1.2675 H   0  0  0  0  0  0  0  0  0  0  0  0
    5.2818   -1.5979    0.3120 H   0  0  0  0  0  0  0  0  0  0  0  0
    4.2343    3.1743   -0.5548 H   0  0  0  0  0  0  0  0  0  0  0  0
    6.8642    0.1634    1.0209 H   0  0  0  0  0  0  0  0  0  0  0  0
    6.3451    2.5508    0.5901 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  8  1  0  0  0  0
  1  9  1  0  0  0  0
  2 12  2  0  0  0  0
  3 15  1  0  0  0  0
  3 34  1  0  0  0  0
  4 15  2  0  0  0  0
  5 16  2  0  0  0  0
  6  8  1  0  0  0  0
  6 10  1  0  0  0  0
  6 12  1  0  0  0  0
  7 11  1  0  0  0  0
  7 16  1  0  0  0  0
  7 33  1  0  0  0  0
  8 11  1  0  0  0  0
  8 24  1  0  0  0  0
  9 10  1  0  0  0  0
  9 13  1  0  0  0  0
  9 14  1  0  0  0  0
 10 15  1  0  0  0  0
 10 25  1  0  0  0  0
 11 12  1  0  0  0  0
 11 26  1  0  0  0  0
 13 27  1  0  0  0  0
 13 28  1  0  0  0  0
 13 29  1  0  0  0  0
 14 30  1  0  0  0  0
 14 31  1  0  0  0  0
 14 32  1  0  0  0  0
 16 17  1  0  0  0  0
 17 18  1  0  0  0  0
 17 35  1  0  0  0  0
 17 36  1  0  0  0  0
 18 19  2  0  0  0  0
 18 20  1  0  0  0  0
 19 21  1  0  0  0  0
 19 37  1  0  0  0  0
 20 22  2  0  0  0  0
 20 38  1  0  0  0  0
 21 23  2  0  0  0  0
 21 39  1  0  0  0  0
 22 23  1  0  0  0  0
 22 40  1  0  0  0  0
 23 41  1  0  0  0  0
M  END
$$$$`;

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
            href="https://github.com/moltools/CineMol/issues" 
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

const CineMol = () => {
    const [isDrawerOpen, setIsDrawerOpen] = useState(false);
    const [statusServer, setStatusServer] = useState(true);

    // General state variables.
    const [version, setVersion] = useState("0.0.0");                    // Version of the CineMol component.
    const [mode, setMode] = useState("dark");                           // Dark or light background of the molecular model.
    const [svgString, setSvgString] = useState("");                     // SVG representation of the molecular model.
    const [isLoading, setIsLoading] = useState(false);                  // App grays out when loading.
    const [initialRender, setInitialRender] = useState(true);           // Initial render of the molecular model.

    // Button state variables.
    const [sdfString, setSdfString] = useState("");                     // SDF string of the molecular model.
    const [style, setStyle] = useState("space-filling");                // Style of the molecular model.
    const [look, setLook] = useState("glossy");                         // Look of the molecular model.
    const [includeHydrogens, setIncludeHydrogens] = useState(false);    // Include hydrogens in the molecular model.
    const [resolution, setResolution] = useState(50);                   // Resolution of the molecular model.
    const [rotationX, setRotationX] = useState(0.0);                    // Rotation over the x-axis of the molecular model in radians.
    const [rotationY, setRotationY] = useState(0.0);                    // Rotation over the y-axis of the molecular model in radians.
    const [rotationZ, setRotationZ] = useState(0.0);                    // Rotation over the z-axis of the molecular model in radians.
    const [viewBox, setViewBox] = useState(undefined);                  // View box of the molecular model.
    const [width, setWidth] = useState(500);                            // Width of the molecular model.
    const [height, setHeight] = useState(500);  
    
    const toggleDrawer = () => {
        setIsDrawerOpen(!isDrawerOpen);
    };

    // Downloads the SVG representation of the molecular model.
    const handleDownloadSvgString = () => {
        if (svgString.length > 0) {
            const blob = new Blob([svgString], { type: "image/svg+xml" });
            const url = window.URL.createObjectURL(blob);
            const a = document.createElement("a");
            a.href = url;
            a.download = "model.svg";
            document.body.appendChild(a);
            a.click();
            window.URL.revokeObjectURL(url);
            document.body.removeChild(a);
        } else {
            toast.error("No SVG model to download!", { autoClose: true });
        }
    };

    // Uploads an SDF file and sets the SDF string.
    const handleUploadSdfFile = () => {
        const fileInput = document.createElement("input");
        fileInput.type = "file";
        fileInput.accept = ".sdf";
        fileInput.onchange = event => {
            const file = event.target.files[0];

            // Check if the file size is greater than 1 MB (1,000,000 bytes).
            if (file.size > 1000000) {
                toast.error("File size exceeds 1 MB limit!", { autoClose: true });
                return;
            };

            const reader = new FileReader();
            reader.onload = event => {
                const sdfString = event.target.result;
                toast.success("SDF file uploaded!", { autoClose: true });
                setSdfString(sdfString);
            };
            reader.readAsText(file);
        };
        fileInput.click();
    };

    // Get the version of the CineMol component.
    const handleGetVersion = async () => {
        try {
            const response = await fetch("/api/fetch_cinemol_version");
            
            if (!response.ok) {
                throw new Error("Network response was not ok!");
            };

            const json = await response.json();

            // Unpack response.
            if (json.status === "success") {
                setVersion(json.payload.version);
            } else if (json.status === "warning") {
                toast.warn(json.message);
            } else if (json.status === "failure") {
                toast.error(json.message);
            };

        } catch (error) {
            const msg = "Could not retrieve version!";
            toast.error(msg, { autoClose: true });
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

    // Draw the molecular model. 
    const handleDrawModel = async () => {
        // Set is loading to true, which grays out model and deactivates buttons.
        setIsLoading(true);

        if (sdfString.length === 0) {
            setIsLoading(false);
            return;
        };
    
        try {
            const response = await fetch("/api/draw_model", {
                method: "POST",
                headers: { "Content-Type": "application/json" },
                body: JSON.stringify({
                    "sdf_string": sdfString,
                    "style": style,
                    "look": look,
                    "include_hydrogens": includeHydrogens,
                    "resolution": resolution,
                    "rotation_x": rotationX,
                    "rotation_y": rotationY,
                    "rotation_z": rotationZ,
                    "view_box": viewBox,
                    "width": width,
                    "height": height,
                }),
            });
    
            if (!response.ok) {
                throw new Error("Network response was not ok!");
            };
    
            const json = await response.json();

            if (json.status === "success") {
                setSvgString(json.payload.svg_string);
                setViewBox(json.payload.view_box);
            } else if (json.status === "warning") {
                toast.warn(json.message);
            } else if (json.status === "failure") {
                toast.error(json.message);
            };
            setIsLoading(false);

        } catch (error) {
            const msg = "Could not parse SDF string!";
            toast.error(msg, { autoClose: true });
            console.error(error);
            setIsLoading(false);
        };
    };

    const handleToggleStyle = () => {
        if        (style === "space-filling")   { setStyle("ball-and-stick");
        } else if (style === "ball-and-stick")  { setStyle("tube");
        } else if (style === "tube")            { setStyle("wireframe");
        } else if (style === "wireframe")       { setStyle("space-filling");
        } else                                  { setStyle("ball-and-stick"); }
    };

    const handleToggleLook = () => {
        if        (look === "glossy")   { setLook("cartoon");
        } else if (look === "cartoon")  { setLook("glossy");
        } else                          { setLook("cartoon"); }
    };

    const handleToggleMode = () => {
        if (mode === "dark") { setMode("light"); } else { setMode("dark"); };
    };

    const handleReset = () => {
        setMode("light");
        setSvgString("");
        setSdfString("");
        setStyle("space-filling");
        setLook("glossy");
        setIncludeHydrogens(false);
        setResolution(50);
        setRotationX(0.0);
        setRotationY(0.0);
        setRotationZ(0.0);
        setViewBox(undefined);
        setWidth(500);
        setHeight(500);
    };

    /**
     * Get the version of the CineMol component and skip the first render. When any of the following state variables 
     * change, the component will re-render and the useEffect hook will run again.
     */ 
    useEffect(() => {
        if (initialRender) {
            handleGetVersion();
            setInitialRender(false);
            // toast.warn("CineMol is currently under development, and may not look as expected!", { autoClose: false });
        } else {
            handleDrawModel();
        };
    }, [sdfString, style, look, includeHydrogens, resolution, rotationX, rotationY, rotationZ, width, height]);

    useEffect(() => {
        if (isDrawerOpen) {
            checkStatusServer();
        };
    }, [isDrawerOpen]);


    return (
        <Box sx={{ display: "flex", width: "100%", flexDirection: "column" }}>
            <CssBaseline />
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
                                <img src="/cinemol.svg" alt="CineMol" style={{ width: 40, height: 40 }} />
                            </Box>
                            <Typography variant="h6" noWrap>
                                {`CineMol (${version})`}
                            </Typography>
                        </Box>
                    </Box>
                </Toolbar>
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
                </List>
            </Drawer>
            </AppBar>
            <Box component="main">
                <Toolbar />
                <Box sx={{
                    m: 0, 
                    p:0, 
                    position: "relative",
                    height: "calc(100vh - 64px)", // 64px is the height of the AppBar.
                    width: "100vw",
                    overflow: "hidden",
                }}> 
                    {isLoading && <LoadingOverlay />}
                    <Box sx={{ 
                        height: "calc(100vh - 64px)", // 64px is the height of the AppBar.
                        width: "100vw",
                        backgroundColor: mode === "dark" ? "#000" : "#f5f5f5", 
                    }}>
                        <div dangerouslySetInnerHTML={{ __html: svgString }} />
                    </Box>
                    <Draggable
                        defaultPosition={{ x: 0, y: 0 }}
                        scale={1}
                        handle=".drag-handle"
                    >
                        <Box sx={{ 
                            position: "absolute",
                            top: 10,
                            left: 10,
                            display: "flex", 
                            flexDirection: "column", 
                            maxWidth: "250px" 
                        }}>
                            <Box 
                                sx={{ 
                                    display: "flex", 
                                    flexDirection: "row", 
                                    justifyContent: "center", 
                                    alignItems: "center", 
                                    backgroundColor: "#555", 
                                    padding: 1, 
                                    borderRadius: 4, 
                                    borderBottomLeftRadius: 0, 
                                    borderBottomRightRadius: 0, 
                                    width: "100%",
                                    cursor: "move"
                                }}
                                className="drag-handle"
                            >
                                <DragHandleIcon />
                            </Box>
                            <Box sx={{
                                maxHeight: "250px",
                                overflow: "auto",
                                backgroundColor: "#f0f0f0",
                                boxShadow: 4,
                                borderRadius: 0,
                            }}>
                                <List sx={{ padding: 0 }}>
                                    <Button 
                                        variant="contained" 
                                        color="primary" 
                                        onClick={handleUploadSdfFile} 
                                        disabled={isLoading}
                                        style={{width: "100%", borderRadius: 0, marginTop: "1px", boxShadow: "none", justifyContent: "flex-start"}}
                                    >
                                        Upload SDF
                                    </Button>
                                    <Button 
                                        variant="contained" 
                                        color="primary" 
                                        onClick={handleDownloadSvgString} 
                                        disabled={isLoading}
                                        style={{width: "100%", borderRadius: 0, marginTop: "1px", boxShadow: "none", justifyContent: "flex-start"}}
                                    >
                                        Download SVG
                                    </Button>
                                    <Button 
                                        variant="contained" 
                                        color="primary" 
                                        onClick={() => setSdfString(exampleSdfString)}
                                        disabled={isLoading}
                                        style={{width: "100%", borderRadius: 0, marginTop: "1px", boxShadow: "none", justifyContent: "flex-start"}}
                                    >
                                        Load example
                                    </Button>
                                    <Button 
                                        variant="contained" 
                                        color="primary" 
                                        onClick={handleReset} 
                                        disabled={isLoading}
                                        style={{width: "100%", borderRadius: 0, marginTop: "1px", boxShadow: "none", justifyContent: "flex-start"}}
                                    >
                                        Reset
                                    </Button>
                                    <Button
                                        variant="contained"
                                        color="primary"
                                        onClick={handleToggleStyle}
                                        disabled={isLoading}
                                        style={{width: "100%", borderRadius: 0, marginTop: "1px", boxShadow: "none", justifyContent: "flex-start", flexDirection: "row"}}
                                    >
                                        <SwitchRight sx={{ paddingRight: 1 }} />
                                        Style: {style}
                                    </Button>
                                    <Button
                                        variant="contained"
                                        color="primary"
                                        onClick={handleToggleLook}
                                        disabled={isLoading}
                                        style={{width: "100%", borderRadius: 0, marginTop: "1px", boxShadow: "none", justifyContent: "flex-start", flexDirection: "row"}}
                                    >
                                        <SwitchRight sx={{ paddingRight: 1 }} />
                                        Look: {look}
                                    </Button>
                                    <Button
                                        variant="contained"
                                        color="primary"
                                        onClick={handleToggleMode}
                                        disabled={isLoading}
                                        style={{width: "100%", borderRadius: 0, marginTop: "1px", boxShadow: "none", justifyContent: "flex-start", flexDirection: "row"}}
                                    >
                                        <SwitchRight sx={{ paddingRight: 1 }} />
                                        Mode: {mode}
                                    </Button>
                                    <Button
                                        variant="contained"
                                        color="primary"
                                        onClick={() => setIncludeHydrogens(!includeHydrogens)}
                                        disabled={isLoading}
                                        style={{width: "100%", borderRadius: 0, marginTop: "1px", boxShadow: "none", justifyContent: "flex-start"}}
                                    >
                                        <SwitchRight sx={{ paddingRight: 1 }} />
                                        Hydrogens: {includeHydrogens ? "Yes" : "No"}
                                    </Button>
                                    <Box sx={{ padding: 0, paddingTop: 1, paddingRight: 2, paddingLeft: 2 }}>
                                        <Typography id="discrete-slider" gutterBottom>
                                            {`Resolution (${resolution})`}
                                        </Typography>
                                        <Slider
                                            aria-labelledby="discrete-slider"
                                            valueLabelDisplay="auto"
                                            step={10}
                                            marks
                                            min={40}
                                            max={100}
                                            value={resolution}
                                            disabled={isLoading}
                                            onChange={(event, value) => setResolution(value)}
                                        />
                                    </Box>
                                    <Divider />
                                    <Box sx={{ padding: 0, paddingTop: 1, paddingRight: 2, paddingLeft: 2 }}>
                                        <Typography id="discrete-slider" gutterBottom>
                                            {`Width (${width} px.)`}
                                        </Typography>
                                        <Slider
                                            aria-labelledby="discrete-slider"
                                            valueLabelDisplay="auto"
                                            step={100}
                                            marks
                                            min={100}
                                            max={1000}
                                            value={width}
                                            disabled={isLoading}
                                            onChange={(event, value) => setWidth(value)}
                                        />
                                    </Box>
                                    <Divider />
                                    <Box sx={{ padding: 0, paddingTop: 1, paddingRight: 2, paddingLeft: 2 }}>
                                        <Typography id="discrete-slider" gutterBottom>
                                            {`Height (${height} px.)`}
                                        </Typography>
                                        <Slider
                                            aria-labelledby="discrete-slider"
                                            valueLabelDisplay="auto"
                                            step={100}
                                            marks
                                            min={100}
                                            max={1000}
                                            value={height}
                                            disabled={isLoading}
                                            onChange={(event, value) => setHeight(value)}
                                        />
                                    </Box>
                                    <Divider />
                                    <Box sx={{ padding: 0, paddingTop: 1, paddingRight: 2, paddingLeft: 2 }}>
                                        <Typography id="discrete-slider" gutterBottom>
                                            {`Rotation X (${rotationX} rad.)`}
                                        </Typography>
                                        <Slider
                                            aria-labelledby="discrete-slider"
                                            valueLabelDisplay="auto"
                                            step={Math.PI / 12}
                                            marks
                                            min={0}
                                            max={2 * Math.PI}
                                            value={rotationX}
                                            disabled={isLoading}
                                            onChange={(event, value) => setRotationX(Math.round(value * 100) / 100)}
                                        />
                                    </Box>
                                    <Divider />
                                    <Box sx={{ padding: 0, paddingTop: 1, paddingRight: 2, paddingLeft: 2 }}>
                                        <Typography id="discrete-slider" gutterBottom>
                                            {`Rotation Y (${rotationY} rad.)`}
                                        </Typography>
                                        <Slider
                                            aria-labelledby="discrete-slider"
                                            valueLabelDisplay="auto"
                                            step={Math.PI / 12}
                                            marks
                                            min={0}
                                            max={2 * Math.PI}
                                            value={rotationY}
                                            disabled={isLoading}
                                            onChange={(event, value) => setRotationY(Math.round(value * 100) / 100)}
                                        />
                                    </Box>
                                    <Divider />
                                    <Box sx={{ padding: 0, paddingTop: 1, paddingRight: 2, paddingLeft: 2 }}>
                                        <Typography id="discrete-slider" gutterBottom>
                                            {`Rotation Z (${rotationZ} rad.)`}
                                        </Typography>
                                        <Slider
                                            aria-labelledby="discrete-slider"
                                            valueLabelDisplay="auto"
                                            step={Math.PI / 12}
                                            marks
                                            min={0}
                                            max={2 * Math.PI}
                                            value={rotationZ}
                                            disabled={isLoading}
                                            onChange={(event, value) => setRotationZ(Math.round(value * 100) / 100)}
                                        />
                                    </Box>
                                </List>
                            </Box>
                            <Box 
                                sx={{ 
                                    display: "flex", 
                                    flexDirection: "row", 
                                    justifyContent: "center", 
                                    alignItems: "center", 
                                    backgroundColor: "#555", 
                                    padding: 1, 
                                    borderRadius: 4, 
                                    borderTopLeftRadius: 0, 
                                    borderTopRightRadius: 0, 
                                    width: "100%",
                                    cursor: "move",
                                    paddingTop: "1px",
                                }}
                                className="drag-handle"
                            >
                                <DragHandleIcon />
                            </Box>
                        </Box>
                    </Draggable>
                </Box>
            </Box>
        </Box>
    );
};

export default CineMol;