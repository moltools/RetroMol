import React, { useRef, useEffect, useState } from "react";
import { Box, Button, Divider, FormControl, FormControlLabel, IconButton, InputLabel, MenuItem, MenuProps, Paper, 
    PaperProps, Select, Switch, TextField, Tooltip, Typography,
    Radiogroup, Radio, RadioGroup
} from "@mui/material";
import { DragDropContext, Droppable, Draggable } from "react-beautiful-dnd";
import TypographyTooltip from "../common/TypographyTooltip";
import RefreshIcon from "@mui/icons-material/Refresh";
import CloseIcon from "@mui/icons-material/Close";
import AddIcon from "@mui/icons-material/Add";
import SearchIcon from "@mui/icons-material/Search";
import { toast } from "react-toastify";
import MultiSelect from "../common/MultiSelect";
import { ExpandMore, ExpandLess } from "@mui/icons-material";

const defaultMotif = {
    motifType: "polyketide",
    polyketideType: "Any",
    polyketideDecor: "Any",
    peptideSource: "Any",
    peptideCid: "Any",
    stereoAlpha: "Any",
    stereoBeta: "Any"
};


const drawSmiles = async (smiles, dimensions) => {
    if (smiles.length === 0) {
        return;
    };

    try {
        const response = await fetch("/api/draw_smiles", {
            method: "POST",
            headers: { "Content-Type": "application/json" },
            body: JSON.stringify({
                "smiles": smiles,
                "dimensions": dimensions
            }),
        });

        if (!response.ok) {
            throw new Error("Network response was not ok!");
        };

        const json = await response.json();

        if (json.status === "success") {
            return json.payload.svgString;
        } else if (json.status === "warning") {
            toast.warn(json.message);
        } else if (json.status === "failure") {
            toast.error(json.message);
        };

    } catch (error) {
        const msg = "Could not draw SMILES string!";
        toast.error(msg, { autoClose: true });
        console.error(error);
    };

    return null;
};


const PolyketideDisplayComponent = ({ motif, dimensions }) => {
    const svgContainerRef = useRef(null);
    const motifStr = `${motif.polyketideType}`

    const restGroup = (motif.polyketideDecor === "Any") ? "*" : `[${motif.polyketideDecor}*]`;

    // string to smiles dict
    const smilesDict = {
        "Any": "*",
        "A": `*SC(=O)C(${restGroup})C(=O)[OH]`,
        "B": `*SC([OH])C(${restGroup})C(=O)[OH]`,
        "C": `*SC=C(${restGroup})C(=O)[OH]`,
        "D": `*SCC(${restGroup})C(=O)[OH]`,
    };

    const smiles = smilesDict[motifStr] || "*";

    useEffect(() => {
        const generateSVG = async () => {
            const svgString = await drawSmiles(smiles, dimensions);
            if (svgString && svgContainerRef.current) {
                svgContainerRef.current.innerHTML = svgString;
            }
        };

        generateSVG();
    }, [smiles, dimensions]);

    return (
        <Box sx={{ height: "100px" }}>
            <div ref={svgContainerRef} />
        </Box>
    );
};

const fetchIsomericSmilesFromPubChem = async (cid) => {
    try {
        const response = await fetch(`https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${cid}/property/IsomericSMILES/JSON`);
        if (!response.ok) {
            throw new Error("Failed to fetch SMILES from PubChem");
        }
        const data = await response.json();
        return data.PropertyTable.Properties[0].IsomericSMILES;
    } catch (error) {
        // console.error("Error fetching SMILES from PubChem:", error);
        return "*";
    }
};

const fetchCanonicalSmilesFromPubChem = async (cid) => {
    try {
        const response = await fetch(`https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${cid}/property/CanonicalSMILES/JSON`);
        if (!response.ok) {
            throw new Error("Failed to fetch SMILES from PubChem");
        }
        const data = await response.json();
        return data.PropertyTable.Properties[0].CanonicalSMILES;
    } catch (error) {
        // console.error("Error fetching SMILES from PubChem:", error);
        return "*";
    }
};

async function fetchTrivialNameFromPubChem(cid) {
    const url = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${cid}/synonyms/JSON`;
    try {
        const response = await fetch(url);
        if (!response.ok) {
            throw new Error('Network response was not ok');
        }
        const data = await response.json();
        return data.InformationList.Information[0].Synonym[0];
    } catch (error) {
        // console.error("Error fetching trivial name:", error);
        return "unavailable";
    }
};

const PubchemDisplayComponent = ({ cid, dimensions }) => {
    const svgContainerRef = useRef(null);
    const [trivialName, setTrivialName] = useState("");

    useEffect(() => {
        const fetchData = async () => {
            let smiles = "*";
            let name = "unavailable";

            if (!isNaN(parseInt(cid)) && parseInt(cid) > 0) {
                smiles = await fetchIsomericSmilesFromPubChem(cid);
                name = await fetchTrivialNameFromPubChem(cid); // Assuming this function exists and fetches the trivial name
            }

            if (smiles === "*") {
                smiles = await fetchCanonicalSmilesFromPubChem(cid);
            }

            console.log(name, smiles);

            const svgString = await drawSmiles(smiles, dimensions);
            if (svgString && svgContainerRef.current) {
                svgContainerRef.current.innerHTML = svgString;
            }

            setTrivialName(name);
        };

        fetchData();
    }, [cid, dimensions]);

    return (
        <Box>
            <Tooltip title={trivialName} arrow>
                <Typography variant="body2" gutterBottom style={{ maxWidth: '170px', overflow: 'hidden', textOverflow: 'ellipsis', whiteSpace: 'nowrap' }}>
                    Name: {trivialName}
                </Typography>
            </Tooltip>
            <Box sx={{ height: "100px" }}>
                <div ref={svgContainerRef} />
            </Box>
        </Box>
    );
};

const QueryForm = ({ results, selectedResultIndex, setSelectedResultIndex, columns, setColumns, submit }) => {
    
    // query settings
    const [queryType, setQueryType] = useState("match"); // match or query
    const [alignmentType, setAlignmentType] = useState("local"); // global or local
    const [gapPenalty, setGapPenalty] = useState(3);
    const [endGapPenalty, setEndGapPenalty] = useState(2);
    const [queryHasLeadingModules, setQueryHasLeadingModules] = useState(false);
    const [queryHasTrailingModules, setQueryHasTrailingModules] = useState(false);

    const [allowStereochemistry, setAllowStereochemistry] = useState(false);
    
    const [allBioactivityLabels, setAllBioactivityLabels] = useState([]);
    const [selectedBioactivityLabels, setSelectedBioactivityLabels] = useState([]);
    const [allOrganismLabels, setAllOrganismLabels] = useState([]);
    const [selectedOrganismLabels, setSelectedOrganismLabels] = useState([]);

    const [maxNumMatches, setMaxNumMatches] = useState(10);
    const [queryAgainstMolecules, setQueryAgainstMolecules] = useState(true);
    const [queryAgainstProtoclusters, setQueryAgainstProtoclusters] = useState(false);
    const [minMatchLength, setMinMatchLength] = useState(1);
    const [maxMatchLength, setMaxMatchLength] = useState(100);

    const [showQueryTypeOptions, setShowQueryTypeOptions] = useState(false);
    const [showQuerySpaceOptions, setShowQuerySpaceOptions] = useState(false);

    const handleRefresh = () => {
        setQueryType("match");
        setAlignmentType("local");
        setAllowStereochemistry(false);
        setGapPenalty(2);
        setEndGapPenalty(1);
        setQueryHasLeadingModules(false);
        setQueryHasTrailingModules(false);
        setSelectedBioactivityLabels([]);
        setSelectedOrganismLabels([]);
        setMaxNumMatches(10);
        setQueryAgainstMolecules(true);
        setQueryAgainstProtoclusters(false);
        setMinMatchLength(1);
        setMaxMatchLength(100);
        setSelectedResultIndex(null);
        
        // when you start querying right away, the selectedResultIndex is null
        if (selectedResultIndex === null) {
            setColumns([]);
            return;
        };
        setColumns(results[selectedResultIndex]);
    };

    const onDragEnd = (result) => {
        const { destination, source } = result;
        if (!destination) {
            return;
        };
        const newColumns = Array.from(columns);
        const [movedItem] = newColumns.splice(source.index, 1);
        newColumns.splice(destination.index, 0, movedItem);
        setColumns(newColumns);
    };

    const handleDeleteColumn = (index) => {
        const newColumns = columns.filter((_, i) => i !== index);
        setColumns(newColumns);
    };


    const handleAddColumn = () => {
        const newColumns = [...columns, [defaultMotif]];
        setColumns(newColumns);
    };

    const addMotifToColumn = (index) => () => {
        const newColumns = Array.from(columns);
        newColumns[index] = [...newColumns[index], defaultMotif];
        setColumns(newColumns);
    };

    const fetchBioactivityLabels = async () => {
        try {
            const response = await fetch("/api/fetch_bioactivity_labels", {
                method: "GET",
                headers: { "Content-Type": "application/json" }
            });

            if (!response.ok) {
                throw new Error("Network response was not ok!");
            };

            const json = await response.json();

            if (json.status === "success") {
                setAllBioactivityLabels(json.payload.bioactivityLabels);
            } else if (json.status === "warning") {
                toast.warn(json.message);
            } else if (json.status === "failure") {
                toast.error(json.message);
            };
        } catch (error) {
            toast.error(error.message);
        };
    };

    const fetchOrganismLabels = async () => {
        try {
            const response = await fetch("/api/fetch_organism_labels", {
                method: "GET",
                headers: { "Content-Type": "application/json" }
            });

            if (!response.ok) {
                throw new Error("Network response was not ok!");
            };

            const json = await response.json();

            if (json.status === "success") {
                setAllOrganismLabels(json.payload.organismLabels);
            } else if (json.status === "warning") {
                toast.warn(json.message);
            } else if (json.status === "failure") {
                toast.error(json.message);
            };
        } catch (error) {
            toast.error(error.message);
        };
    };

    useEffect(() => {
        fetchBioactivityLabels();
        fetchOrganismLabels();
    }, []);

    useEffect(() => {
        if (selectedResultIndex !== null) {
            const selectedResult = results[selectedResultIndex]["query"];
            setColumns(selectedResult);
        } else {
            setColumns([]);
        };
    }, [selectedResultIndex]);

    return (
        <Paper
            sx={{
                p: 3,
                backgroundColor: "#f0f0f0",
                borderRadius: 2,
                boxShadow: 3,
                width: "100%",
                border: "1px solid #ccc",
                display: "flex",
                flexDirection: "column",
            }}
        >
            <Box sx={{ flexShrink: 0 }}>
                <Box sx={{ display: "flex", justifyContent: "space-between", alignItems: "center" }}>
                    <Box sx={{ display: "flex", gap: 1 }}>
                        <TypographyTooltip title="The query is used to search for similar protoclusters through pairwise matching or pattern matching.">
                            <Typography variant="h6" gutterBottom>
                                Build query
                            </Typography>
                        </TypographyTooltip>
                    </Box>
                    <Box>
                        <IconButton
                            type="button"
                            onClick={handleRefresh}
                            color="primary"
                            sx={{ transform: "translateY(-3px)" }}
                        >
                            <RefreshIcon />
                        </IconButton>
                    </Box>
                </Box>
                <Divider sx={{ mb: 2 }} />
                <Box>
                    <Box 
                        style={{
                            minHeight: "50px", 
                            backgroundColor: "#555",
                            borderTopLeftRadius: "4px",
                            borderTopRightRadius: "4px",
                        }}
                    >
                        <Button
                            variant="contained"
                            color="primary"
                            sx={{ margin: "8px", width: "150px" }}
                            onClick={handleAddColumn}
                        >
                            Add motif
                        </Button>
                        <Tooltip title="Under development." arrow>
                            <FormControlLabel
                                control={<Switch 
                                    checked={allowStereochemistry}
                                    onChange={(event) => setAllowStereochemistry(event.target.checked)}
                                    color="primary"
                                    name="allowStereochemistry"
                                    disabled={true}
                                />}
                                label={allowStereochemistry ? "PK stereochemistry" : "No PK stereochemistry"}
                                style={{ color: "white", margin: "8px" }}
                            />
                        </Tooltip>
                    </Box>
                    <Box sx={{ 
                        flexGrow: 1, 
                        overflowY: "auto", 
                        flex: "1 1 50%",
                        minHeight: "275px",
                        backgroundColor: "white",
                        boxShadow: "inset 0 2px 5px rgba(0, 0, 0, 0.1)",
                        border: "1px solid #bbb", 
                        borderTopLeftRadius: "0",
                        borderTopRightRadius: "0",
                        borderBottomLeftRadius: "4px",
                        borderBottomRightRadius: "4px",
                        pt: 2,
                    }}>
                        {(columns && columns.length) ? (
                            <Box 
                                sx={{ 
                                    overflowX: "auto", 
                                    minHeight:"255px", 
                                    width: "100%", 
                                    pb: 3, 
                                }}
                            >
                                <DragDropContext onDragEnd={onDragEnd} style={{ height: "100%" }}>
                                    <Droppable droppableId="droppable" direction="horizontal">
                                        {(provided) => (
                                        <div ref={provided.innerRef} {...provided.droppableProps} style={{ display: "flex", width: "max-content", height: "100%" }}>
                                            {columns.map((column, index) => (
                                            <Draggable 
                                                key={index} 
                                                draggableId={`column-${index}`} 
                                                index={index}
                                            >
                                                {(provided) => (
                                                <div 
                                                    ref={provided.innerRef} 
                                                    {...provided.draggableProps} 
                                                    {...provided.dragHandleProps}
                                                    key={`column-${index}`}
                                                >
                                                    <Paper 
                                                        elevation={3} 
                                                        sx={{ backgroundColor: "#f0f0f0" }}
                                                        style={{ 
                                                            minWidth: 200, 
                                                            maxWidth: 200,
                                                            margin: "0 8px",
                                                            marginRight: index === columns.length - 1 ? "8px" : "4px" 
                                                        }}
                                                    >   
                                                        <Box
                                                            sx={{
                                                                backgroundColor: "#555",
                                                                color: "#fff",
                                                                padding: "4px 8px",
                                                                display: "flex",
                                                                justifyContent: "space-between",
                                                                alignItems: "center",
                                                                marginBottom: 1,
                                                                borderTopLeftRadius: "4px",
                                                                borderTopRightRadius: "4px",
                                                                borderBottom: "1px solid #ccc",
                                                            }}
                                                        >
                                                            <Typography variant="h6" component="div">
                                                                Motif {index + 1}
                                                            </Typography>
                                                            <Box sx={{ display: "flex", gap: "3px" }}>
                                                                <Tooltip title="Delete motif.">
                                                                    <IconButton 
                                                                        onClick={() => handleDeleteColumn(index)}
                                                                        style={{ 
                                                                            color: 
                                                                            "#990000", 
                                                                            backgroundColor: "#ff5f57", 
                                                                            padding: "1px",
                                                                            transform: "translateY(-5px) translateX(2px)"
                                                                        }}
                                                                    >
                                                                        <CloseIcon style={{fontSize: "12px"}} />
                                                                    </IconButton>
                                                                </Tooltip>
                                                            </Box>
                                                        </Box>
                                                        <Box style={{ width: "100%" }}>
                                                            {column.map((motif, columnItemIndex) => (
                                                                <Box 
                                                                    sx={{backgroundColor: "#ccc", marginBottom: "10px", minHeight: "130px"}}
                                                                    key={`motif-${selectedResultIndex}-${index}-${columnItemIndex}`}
                                                                >
                                                                    <Box 
                                                                        sx={{
                                                                            backgroundColor: "#ccc",
                                                                            display: "flex",
                                                                            justifyContent: "space-between",
                                                                            alignItems: "center",
                                                                            marginBottom: 1,
                                                                            padding: "4px 8px",
                                                                            paddingTop: "0px",
                                                                            width: "100%",
                                                                        }}
                                                                    >
                                                                        <FormControlLabel
                                                                            control={<Switch
                                                                                checked={motif.motifType === "polyketide"}
                                                                                onChange={(event) => {
                                                                                    const newColumns = columns.map((col, colIndex) => 
                                                                                        colIndex === index 
                                                                                            ? col.map((item, itemIndex) => 
                                                                                                itemIndex === columnItemIndex 
                                                                                                    ? { ...item, motifType: event.target.checked ? "polyketide" : "peptide" } 
                                                                                                    : item
                                                                                            )
                                                                                            : col
                                                                                    );
                                                                                    setColumns(newColumns);
                                                                                }}
                                                                                color="primary"
                                                                                name={`motif-${index}-${columnItemIndex}`}
                                                                            />}
                                                                            label={motif.motifType === "peptide" ? "Other" : motif.motifType.charAt(0).toUpperCase() + motif.motifType.slice(1)}
                                                                        />
                                                                        <Box sx={{ display: "flex", gap: "3px" }}>
                                                                            {(motif.motifType === "peptide") && (
                                                                                <Tooltip title="Search for ID on PubChem.">
                                                                                    <IconButton 
                                                                                        onClick={() => {
                                                                                            // check if the PubChem ID is valid by seeing if you can parse it as integer
                                                                                            if (isNaN(parseInt(motif.peptideCid))) {
                                                                                                toast.warn("Invalid PubChem ID.");
                                                                                                return;
                                                                                            };
                                                                                            window.open(`https://pubchem.ncbi.nlm.nih.gov/compound/${motif.peptideCid}`, "_blank");
                                                                                        }}
                                                                                        style={{ 
                                                                                            color: "#985600", 
                                                                                            backgroundColor: "#febc2e", 
                                                                                            padding: "1px",
                                                                                            transform: "translateY(-5px) translateX(2px)"
                                                                                        }}
                                                                                    >
                                                                                        <SearchIcon style={{fontSize: "12px" }} />
                                                                                    </IconButton>
                                                                                </Tooltip>
                                                                            )}
                                                                            {(column.length > 1) && (
                                                                                <Tooltip title="Delete option.">
                                                                                    <IconButton 
                                                                                        onClick={() => {
                                                                                            const newColumns = columns.map((col, colIndex) => 
                                                                                                colIndex === index 
                                                                                                    ? col.filter((_, itemIndex) => itemIndex !== columnItemIndex)
                                                                                                    : col
                                                                                            );
                                                                                            setColumns(newColumns);
                                                                                        }}
                                                                                        style={{ 
                                                                                            color: "#990000", 
                                                                                            backgroundColor: "#ff5f57", 
                                                                                            padding: "1px",
                                                                                            transform: "translateY(-5px) translateX(2px)"
                                                                                        }}
                                                                                    >
                                                                                        <CloseIcon style={{fontSize: "12px"}} />
                                                                                    </IconButton>
                                                                                </Tooltip>
                                                                            )}
                                                                        </Box>
                                                                    </Box>
                                                                    {motif.motifType === "polyketide" && (
                                                                        <Box sx={{ display: "flex", flexDirection: "column", alignItems: "center", justifyContent: "space-between", gap: "20px", paddingTop: "25px" }}>
                                                                            <PolyketideDisplayComponent
                                                                                motif={motif}
                                                                                dimensions={{ width: 170, height: 100 }}
                                                                                key={`polyketide-${selectedResultIndex}-${index}-${columnItemIndex}`}
                                                                            />
                                                                            <Box sx={{ display: "flex", flexDirection: "row", alignItems: "center", justifyContent: "space-between", paddingBottom: "20px" }}>
                                                                                <Box sx={{ width: "100%", paddingLeft: "10px", paddingRight: "2.5px" }}>
                                                                                    <FormControl
                                                                                        variant="filled"
                                                                                        style={{ 
                                                                                            width: "100%",
                                                                                            borderRadius: "0px",
                                                                                        }}
                                                                                    >
                                                                                        <InputLabel id="demo-simple-select-label">Type</InputLabel>
                                                                                        <Select
                                                                                            labelId="demo-simple-select-label"
                                                                                            id="demo-simple-select"
                                                                                            value={motif.polyketideType}
                                                                                            onChange={(event) => {
                                                                                                const newColumns = columns.map((col, colIndex) => 
                                                                                                    colIndex === index 
                                                                                                        ? col.map((item, itemIndex) => 
                                                                                                            itemIndex === columnItemIndex 
                                                                                                                ? { ...item, polyketideType: event.target.value } 
                                                                                                                : item
                                                                                                        )
                                                                                                        : col
                                                                                                );
                                                                                                setColumns(newColumns);
                                                                                            }}
                                                                                            style={{ width: "80px", padding: "4px" }}
                                                                                            MenuProps={{
                                                                                                PaperProps: {
                                                                                                    style: {
                                                                                                        maxHeight: 300,
                                                                                                        overflowY: "auto",
                                                                                                    },
                                                                                                },   
                                                                                            }}
                                                                                        >
                                                                                            <MenuItem value="Any">Any</MenuItem>
                                                                                            <MenuItem value="A">A</MenuItem>
                                                                                            <MenuItem value="B">B</MenuItem>
                                                                                            {allowStereochemistry && <MenuItem value="Br">Br</MenuItem>}
                                                                                            {allowStereochemistry && <MenuItem value="Bs">Bs</MenuItem>}
                                                                                            <MenuItem value="C">C</MenuItem>
                                                                                            <MenuItem value="D">D</MenuItem>
                                                                                            {allowStereochemistry && <MenuItem value="Dr">Dr</MenuItem>}
                                                                                            {allowStereochemistry && <MenuItem value="Ds">Ds</MenuItem>}
                                                                                        </Select>
                                                                                    </FormControl>
                                                                                </Box>
                                                                                <Box sx={{ width: "100%", paddingRight: "10px", paddingLeft: "2.5px" }}>
                                                                                    <FormControl
                                                                                        variant="filled"
                                                                                        style={{
                                                                                            width: "100%",
                                                                                            borderRadius: "0px",
                                                                                        }}
                                                                                    >
                                                                                        <InputLabel id="demo-simple-select-label">Decor</InputLabel>
                                                                                        <Select
                                                                                            labelId="demo-simple-select-label"
                                                                                            id="demo-simple-select"
                                                                                            value={motif.polyketideDecor}
                                                                                            onChange={(event) => {
                                                                                                const newColumns = columns.map((col, colIndex) =>   
                                                                                                    colIndex === index
                                                                                                        ? col.map((item, itemIndex) =>
                                                                                                            itemIndex === columnItemIndex
                                                                                                                ? { ...item, polyketideDecor: event.target.value }
                                                                                                                : item
                                                                                                        )
                                                                                                        : col
                                                                                                );
                                                                                                setColumns(newColumns);
                                                                                            }}
                                                                                            style={{ width: "80px", padding: "4px" }}
                                                                                            MenuProps={{
                                                                                                PaperProps: {
                                                                                                    style: {
                                                                                                        maxHeight: 300,
                                                                                                        overflowY: "auto",
                                                                                                    },
                                                                                                },   
                                                                                            }}
                                                                                        >
                                                                                            <MenuItem value="Any">Any</MenuItem>
                                                                                            <MenuItem value="1">1</MenuItem>
                                                                                            {allowStereochemistry && <MenuItem value="1s">1s</MenuItem>}
                                                                                            {allowStereochemistry && <MenuItem value="1r">1r</MenuItem>}
                                                                                            <MenuItem value="2">2</MenuItem>
                                                                                            {allowStereochemistry && <MenuItem value="2s">2s</MenuItem>}
                                                                                            {allowStereochemistry && <MenuItem value="2r">2r</MenuItem>}
                                                                                            <MenuItem value="3">3</MenuItem>
                                                                                            {allowStereochemistry && <MenuItem value="3s">3s</MenuItem>}
                                                                                            {allowStereochemistry && <MenuItem value="3r">3r</MenuItem>}
                                                                                            <MenuItem value="4">4</MenuItem>
                                                                                            {allowStereochemistry && <MenuItem value="4s">4s</MenuItem>}
                                                                                            {allowStereochemistry && <MenuItem value="4r">4r</MenuItem>}
                                                                                            <MenuItem value="5">5</MenuItem>
                                                                                            {allowStereochemistry && <MenuItem value="5s">5s</MenuItem>}
                                                                                            {allowStereochemistry && <MenuItem value="5r">5r</MenuItem>}
                                                                                            <MenuItem value="6">6</MenuItem>
                                                                                            {allowStereochemistry && <MenuItem value="6s">6s</MenuItem>}
                                                                                            {allowStereochemistry && <MenuItem value="6r">6r</MenuItem>}
                                                                                            <MenuItem value="7">7</MenuItem>
                                                                                            {allowStereochemistry && <MenuItem value="7s">7s</MenuItem>}
                                                                                            {allowStereochemistry && <MenuItem value="7r">7r</MenuItem>}
                                                                                            <MenuItem value="8">8</MenuItem>
                                                                                            {allowStereochemistry && <MenuItem value="8s">8s</MenuItem>}
                                                                                            {allowStereochemistry && <MenuItem value="8r">8r</MenuItem>}
                                                                                            <MenuItem value="9">9</MenuItem>
                                                                                            {allowStereochemistry && <MenuItem value="9s">9s</MenuItem>}
                                                                                            {allowStereochemistry && <MenuItem value="9r">9r</MenuItem>}
                                                                                            <MenuItem value="10">10</MenuItem>
                                                                                            {allowStereochemistry && <MenuItem value="10s">10s</MenuItem>}
                                                                                            {allowStereochemistry && <MenuItem value="10r">10r</MenuItem>}
                                                                                            <MenuItem value="11">11</MenuItem>
                                                                                            {allowStereochemistry && <MenuItem value="11s">11s</MenuItem>}
                                                                                            {allowStereochemistry && <MenuItem value="11r">11r</MenuItem>}
                                                                                            <MenuItem value="12">12</MenuItem>
                                                                                            {allowStereochemistry && <MenuItem value="12s">12s</MenuItem>}
                                                                                            {allowStereochemistry && <MenuItem value="12r">12r</MenuItem>}
                                                                                            <MenuItem value="13">13</MenuItem>
                                                                                            {allowStereochemistry && <MenuItem value="13s">13s</MenuItem>}
                                                                                            {allowStereochemistry && <MenuItem value="13r">13r</MenuItem>}
                                                                                            <MenuItem value="14">14</MenuItem>
                                                                                            {allowStereochemistry && <MenuItem value="14s">14s</MenuItem>}
                                                                                            {allowStereochemistry && <MenuItem value="14r">14r</MenuItem>}
                                                                                            <MenuItem value="15">15</MenuItem>
                                                                                            {allowStereochemistry && <MenuItem value="15s">15s</MenuItem>}
                                                                                            {allowStereochemistry && <MenuItem value="15r">15r</MenuItem>}
                                                                                        </Select>
                                                                                    </FormControl>
                                                                                </Box>
                                                                            </Box>
                                                                        </Box>
                                                                    )}
                                                                    {motif.motifType === "peptide" && (
                                                                        <Box sx={{ display: "flex", flexDirection: "column", alignItems: "center", justifyContent: "space-between", gap: "20px" }}>
                                                                            <PubchemDisplayComponent
                                                                                cid={motif.peptideCid}
                                                                                dimensions={{ width: 170, height: 100 }}
                                                                                key={`peptide-${selectedResultIndex}-${index}-${columnItemIndex}`}
                                                                            />
                                                                            <Box sx={{ display: "flex", flexDirection: "row", alignItems: "center", justifyContent: "space-between", paddingBottom: "20px"}}>
                                                                                <Box sx={{ width: "100%", paddingRight: "10px", paddingLeft: "10px" }}>
                                                                                    <TextField
                                                                                        value={motif.peptideCid} 
                                                                                        defaultValue={"Any"}
                                                                                        onChange={(event) => {
                                                                                            const newColumns = columns.map((col, colIndex) =>
                                                                                                colIndex === index
                                                                                                    ? col.map((item, itemIndex) =>
                                                                                                        itemIndex === columnItemIndex
                                                                                                            ? { ...item, peptideCid: event.target.value }
                                                                                                            : item
                                                                                                    )
                                                                                                    : col
                                                                                            );
                                                                                            setColumns(newColumns);
                                                                                        }}
                                                                                        label="PubChem ID"
                                                                                        style={{ width: "100%", padding: "4px" }}
                                                                                    />
                                                                                </Box>
                                                                            </Box>
                                                                        </Box>
                                                                    )}
                                                                </Box>
                                                            ))}
                                                        </Box>
                                                        <Tooltip title="Add option">
                                                            <Box
                                                                type="button"
                                                                sx={{
                                                                    backgroundColor: "#1976d2",
                                                                    color: "white",
                                                                    height: "35px",
                                                                    display: "flex",
                                                                    alignItems: "center",
                                                                    justifyContent: "center",
                                                                    marginBottom: 1,
                                                                    borderBottomLeftRadius: "4px",
                                                                    borderBottomRightRadius: "4px",
                                                                    cursor: "pointer",
                                                                    width: "100%",
                                                                }}
                                                                onClick={addMotifToColumn(index)}
                                                            >
                                                                <AddIcon style={{ alignSelf: "center" }} />
                                                            </Box>
                                                        </Tooltip>
                                                    </Paper>
                                                </div>
                                                )}
                                            </Draggable>
                                            ))}
                                            {provided.placeholder}
                                        </div>
                                        )}
                                    </Droppable>
                                </DragDropContext>
                            </Box>
                        ) : (
                            <Typography variant="body1" sx={{ mb: 2, paddingLeft: "10px" }}>
                                {`Query is empty.`}
                            </Typography>
                        )}
                    </Box>
                </Box>

                <Box sx={{ display: "flex", justifyContent: "space-between", alignItems: "center", marginTop: "20px" }}>
                    <Box sx={{ display: "flex", gap: 1 }}>
                        <TypographyTooltip title="Choose between pairwise matching or exact pattern matching.">
                            <Typography variant="h6" gutterBottom>
                                Query type
                            </Typography>
                        </TypographyTooltip>
                    </Box>
                </Box>
                <Divider sx={{ mb: 2 }} />
                <Box sx={{ 
                    flexGrow: 1, 
                    flex: "1 1 50%",
                    paddingLeft: "10px"
                }}>
                    <RadioGroup
                        value={queryType}
                        onChange={(event) => setQueryType(event.target.value)}
                        sx={{ mb: 2 }}
                    >
                        <FormControlLabel
                            value="match"
                            control={<Radio />}
                            label="Pairwise matching"
                        />
                        <FormControlLabel
                            value="query"
                            control={<Radio />}
                            label="Pattern matching"
                        />
                    </RadioGroup>
                </Box>
                <Box sx={{ 
                    display: "flex", 
                    flexDirection: "column", 
                    border: "none",
                    borderRadius: "4px" 
                }}>
                {/* <Box> */}
                    <Button
                        variant="outlined"
                        color="primary"
                        onClick={() => setShowQueryTypeOptions(!showQueryTypeOptions)}
                        style={{ 
                            width: "100%", 
                            display: "flex", 
                            justifyContent: "space-between",
                            borderBottomLeftRadius: showQueryTypeOptions ? "0px" : "4px",
                            borderBottomRightRadius: showQueryTypeOptions ? "0px" : "4px",
                        }}
                        sx={{
                            borderColor: "primary.main",
                        
                        }}
                    >
                        {showQueryTypeOptions ? <ExpandLess /> : <ExpandMore />}
                        {showQueryTypeOptions ? "Hide options" : "Show options"}
                        {showQueryTypeOptions ? <ExpandLess /> : <ExpandMore />}
                    </Button>
                {/* </Box> */}
                {showQueryTypeOptions && (
                <Box sx={{ 
                    flexGrow: 1, 
                    // flex: "1 1 50%",
                    padding: "10px",
                    border: "1px solid #ccc",
                    borderColor: "primary.main",
                    backgroundColor: "#cccccc",
                    borderTop: "none",
                    borderBottomLeftRadius: "4px",
                    borderBottomRightRadius: "4px",
                }}>
                    {queryType === "match" && (
                        <RadioGroup
                            value={alignmentType}
                            onChange={(event) => setAlignmentType(event.target.value)}
                            sx={{ mb: 2 }}
                        >
                            <FormControlLabel
                                value="local"
                                control={<Radio />}
                                label="Local alignment strategy"
                            />
                            <FormControlLabel
                                value="global"
                                control={<Radio />}
                                label="Global alignment strategy"
                            />
                        </RadioGroup>
                    )}
                    {queryType === "match" && (
                        <Box>
                            <TextField
                                value={gapPenalty}
                                onChange={(event) => {
                                    const value = event.target.value;
                                    const intValue = parseInt(value);
                                    if (isNaN(intValue) || intValue < 1) {
                                        return;
                                    };
                                    setGapPenalty(intValue);
                                }}
                                label="Gap penalty"
                                margin="normal"
                                variant="outlined"
                                type="number"
                                sx={{ width: "100%", mb: 2 }}
                            />
                            <TextField
                                value={endGapPenalty}
                                onChange={(event) => {
                                    const value = event.target.value;
                                    const intValue = parseInt(value);
                                    if (isNaN(intValue) || intValue < 1) {
                                        return;
                                    };
                                    setEndGapPenalty(intValue);
                                }}
                                label="End gap penalty"
                                margin="normal"
                                variant="outlined"
                                type="number"
                                sx={{ width: "100%" }}
                            />
                        </Box>
                    )}
                    {queryType === "query" && (
                        <Box sx={{ display: "flex", flexDirection: "column", gap: 2 }}>
                            <FormControlLabel
                                control={<Switch 
                                    checked={queryHasLeadingModules}
                                    onChange={(event) => setQueryHasLeadingModules(event.target.checked)}
                                    color="primary"
                                    name="queryHasLeadingModules"
                                />}
                                label="Query has leading modules"
                            />
                            <FormControlLabel
                                control={<Switch 
                                    checked={queryHasTrailingModules}
                                    onChange={(event) => setQueryHasTrailingModules(event.target.checked)}
                                    color="primary"
                                    name="queryHasTrailingModules"
                                />}
                                label="Query has trailing modules"
                            />
                        </Box>
                    )}
                </Box>
                )}
                </Box>

                <Box sx={{ display: "flex", justifyContent: "space-between", alignItems: "center", marginTop: "20px" }}>
                    <Box sx={{ display: "flex", gap: 1 }}>
                        <TypographyTooltip title="Select filters to filter results on. If no filters are selected, all results will be shown. Setting bioactivity or genus filters will filter out all protoclusters.">
                            <Typography variant="h6" gutterBottom>
                                Filters
                            </Typography>
                        </TypographyTooltip>
                    </Box>
                </Box>
                <Divider sx={{ mb: 4 }} />
                {/* <Box 
                    sx={{ 
                        display: "flex",
                        flexDireciton: "row",
                        gap: 2,
                        padding: "10px"
                    }}
                > */}
                    <Box sx={{ m: 0, p: 0, mb: 2, marginLeft: "-5px", paddingRight: "10px" }}>
                    <MultiSelect 
                        title="Select for bioactivity"
                        labels={allBioactivityLabels}
                        selectedLabels={selectedBioactivityLabels}
                        setSelectedLabels={setSelectedBioactivityLabels}
                    />
                    </Box>
                    <Box sx={{ m: 0, p: 0, mb: 2, marginLeft: "-5px", paddingRight: "10px" }}>
                    <MultiSelect 
                        title="Select for genus"
                        labels={allOrganismLabels}
                        selectedLabels={selectedOrganismLabels}
                        setSelectedLabels={setSelectedOrganismLabels}
                    />
                    </Box>
                {/* </Box> */}

                <Box sx={{ display: "flex", justifyContent: "space-between", alignItems: "center", marginTop: "20px" }}>
                    <Box sx={{ display: "flex", gap: 1 }}>
                        <TypographyTooltip title="Further define the scope of the query.">
                            <Typography variant="h6" gutterBottom>
                                Query space
                            </Typography>
                        </TypographyTooltip>
                    </Box>
                </Box>
                <Divider sx={{ mb: 2 }} />
                <Box sx={{ 
                    flexGrow: 1, 
                    flex: "1 1 50%",
                    padding: "10px"
                }}>
                    <Box sx={{ display: "flex", flexDirection: "column", gap: 2 }}>
                        <FormControlLabel
                            control={<Switch 
                                checked={queryAgainstMolecules}
                                onChange={(event) => setQueryAgainstMolecules(event.target.checked)}
                                color="primary"
                                name="queryAgainstMolecules"
                            />}
                            label="Query against molecules"
                        />
                        <FormControlLabel
                            control={<Switch 
                                checked={queryAgainstProtoclusters}
                                onChange={(event) => setQueryAgainstProtoclusters(event.target.checked)}
                                color="primary"
                                name="queryAgainstProtoclusters"
                            />}
                            label="Query against protoclusters"
                        />
                    </Box>
                </Box>
                <Box sx={{ 
                    display: "flex", 
                    flexDirection: "column", 
                    border: "none",
                    borderRadius: "4px" ,
                    paddingTop: "5px",
                }}>
                    <Button
                            variant="outlined"
                            color="primary"
                            onClick={() => setShowQuerySpaceOptions(!showQuerySpaceOptions)}
                            style={{ 
                                width: "100%", 
                                display: "flex", 
                                justifyContent: "space-between",
                                borderBottomLeftRadius: showQuerySpaceOptions ? "0px" : "4px",
                                borderBottomRightRadius: showQuerySpaceOptions ? "0px" : "4px",
                            }}
                            sx={{
                                borderColor: "primary.main",
                            
                            }}
                        >
                            {showQuerySpaceOptions ? <ExpandLess /> : <ExpandMore />}
                            {showQuerySpaceOptions ? "Hide options" : "Show options"}
                            {showQuerySpaceOptions ? <ExpandLess /> : <ExpandMore />}
                    </Button>
                    {showQuerySpaceOptions && (
                    // <Box sx={{ 
                    //     flexGrow: 1, 
                    //     flex: "1 1 50%",
                    //     padding: "10px"
                    // }}>
                    <Box sx={{ 
                        flexGrow: 1, 
                        // flex: "1 1 50%",
                        padding: "10px",
                        border: "1px solid #ccc",
                        borderColor: "primary.main",
                        backgroundColor: "#cccccc",
                        borderTop: "none",
                        borderBottomLeftRadius: "4px",
                        borderBottomRightRadius: "4px",
                    }}>
                        <TextField
                            value={maxNumMatches}
                            onChange={(event) => {
                                const value = event.target.value;
                                const intValue = parseInt(value);
                                if (isNaN(intValue) || intValue < 1 || intValue > 100) {
                                    return;
                                };
                                setMaxNumMatches(intValue);
                            
                            }}
                            label="Max shown matches"
                            margin="normal"
                            variant="outlined"
                            type="number"
                            sx={{ width: "100%", mb: 2 }}
                        />
                        <TextField
                            value={minMatchLength}
                            onChange={(event) => {
                                const value = event.target.value;
                                const intValue = parseInt(value);
                                if (isNaN(intValue) || intValue < 1 || intValue > 100 || intValue > maxMatchLength) {
                                    return;
                                };
                                setMinMatchLength(intValue);
                            }}
                            label="Min length motif code"
                            margin="normal"
                            variant="outlined"
                            type="number"
                            sx={{ width: "100%", mb: 2 }}
                        />
                        <TextField
                            value={maxMatchLength}
                            onChange={(event) => {
                                const value = event.target.value;
                                const intValue = parseInt(value);
                                if (isNaN(intValue) || intValue < 1 || intValue > 100 || intValue < minMatchLength) {
                                    return;
                                };
                                setMaxMatchLength(intValue);
                            }}
                            label="Max length motif code"
                            margin="normal"
                            variant="outlined"
                            type="number"
                            sx={{ width: "100%", mb: 2 }}
                        />
                    </Box>
                    )}
                </Box>
                
                <Box sx={{ display: "flex", justifyContent: "center", alignItems: "center", marginTop: "20px" }}>
                    <Button
                        variant="contained"
                        color="primary"
                        onClick={() => 
                            submit({
                                query: columns,
                                queryType: queryType,
                                alignmentType: alignmentType,
                                gapPenalty: gapPenalty,
                                endGapPenalty: endGapPenalty,
                                queryHasLeadingModules: queryHasLeadingModules,
                                queryHasTrailingModules: queryHasTrailingModules,
                                selectedBioactivityLabels: selectedBioactivityLabels,
                                selectedOrganismLabels: selectedOrganismLabels,
                                maxNumMatches: maxNumMatches,
                                queryAgainstMolecules: queryAgainstMolecules,
                                queryAgainstProtoclusters: queryAgainstProtoclusters,
                                minMatchLength: minMatchLength,
                                maxMatchLength: maxMatchLength,
                            })
                        }
                        style={{ width: "100%" }}
                    >
                        Submit query
                    </Button>
                </Box>

            </Box>
        </Paper>
    );
};

export default QueryForm;