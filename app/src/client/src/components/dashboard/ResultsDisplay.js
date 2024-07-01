import React, { useRef, useEffect, useState } from "react";
import {
    Box,
    Button,
    Divider,
    IconButton,
    Paper,
    Typography
} from "@mui/material";
import TypographyTooltip from "../common/TypographyTooltip";
import RefreshIcon from "@mui/icons-material/Refresh";
import { toast } from "react-toastify";

const ResultsDisplay = ({ 
    results, 
    setResults,
    selectedResultIndex,
    setSelectedResultIndex
}) => {
    const boxRef = useRef(null);
    const [dimensions, setDimensions] = useState({ width: 0, height: 522 });
    const [svgString, setSvgString] = useState("");
    const [smiles, setSmiles] = useState("");
    const [highlights, setHighlights] = useState([]);

    useEffect(() => {
        function updateSize() {
            if (boxRef.current) {
                setDimensions({
                    width: boxRef.current.offsetWidth,
                    // height: boxRef.current.offsetHeight,
                    height: 522
                });
            }
        }
    
        window.addEventListener('resize', updateSize);
        updateSize(); // Initial size update
    
        return () => window.removeEventListener('resize', updateSize);
    }, []); // Empty array ensures this effect runs only once at mount
    

    const handleRefresh = () => {
        setSelectedResultIndex(null);
        setResults([]);
    };

    const handleDrawSmilesWithHighlights = async () => { 

        // also check if it is either null or undefined
        if (smiles === null || smiles === undefined) {
            return;
        }

        // check if smiles is defined
        if (smiles === "") {
            return;
        }

        if (smiles.length === 0) {
            return;
        };
    
        try {
            // console.log(dimensions)
            const response = await fetch("/api/draw_smiles_with_highlights", {
                method: "POST",
                headers: { "Content-Type": "application/json" },
                body: JSON.stringify({
                    "smiles": smiles,
                    "highlights": highlights,
                    "dimensions": dimensions
                }),
            });
    
            if (!response.ok) {
                throw new Error("Network response was not ok!");
            };
    
            const json = await response.json();

            if (json.status === "success") {
                setSvgString(json.payload.svgString);
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
    };

    useEffect (() => {

        // if (prevDimensions.width === dimensions.width && prevDimensions.height === dimensions.height) {
        //     return;
        // };

        // console.log(prevDimensions.width === dimensions.width && prevDimensions.height === dimensions.height)
        // console.log(prevDimensions, dimensions);
        // console.log("Drawing smiles with highlights...")
        
        handleDrawSmilesWithHighlights();
    }, [smiles, highlights, dimensions]);

    const handleSelectResultIndex = (index) => {
        if (selectedResultIndex === index) {
            setSelectedResultIndex(null); // Deselect the result if it is already selected.
        } else {
            setSelectedResultIndex(index); // Select the result.
            setSmiles(results[index].metaData.inputSmiles);
            setHighlights(results[index].metaData.mapping)
        }
    };

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
                // minHeight: "385px",
                // maxHeight: "400px",
            }}
        >
            <Box sx={{ flexShrink: 0 }}>
                <Box sx={{ display: "flex", justifyContent: "space-between", alignItems: "center" }}>
                    <Box sx={{ display: "flex", gap: 1 }}>
                        <TypographyTooltip title="One result resembles a single putative protocluster.">
                            <Typography variant="h6" gutterBottom>
                                Select result
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
            </Box>
            <Box sx={{ display: "flex", flexDirection: "column", gap: 2 }}>
                <Box 
                    ref={boxRef}
                    sx={{ 
                        flex: "1 1 100%",
                        display: "flex", 
                        alignItems: "left", 
                        justifyContent: "left", 
                        padding: "10px", // delete this line when centering component to be implemented
                        backgroundColor: "white", 
                        borderRadius: "4px", 
                        border: "1px solid #bbb",
                        boxShadow: "inset 0 2px 5px rgba(0, 0, 0, 0.1)",
                    }}
                >
                    {selectedResultIndex !== null ? (
                        <Box>
                            {results[selectedResultIndex].queryType === "retrosynthesis" ? (
                                <Box sx={{
                                    height: "500px",
                                    display: "flex", 
                                    alignItems: "center", 
                                    justifyContent: "center", 
                                }}>
                                    {svgString !== "" ? (
                                        <Box sx={{ 
                                            height: "500px",
                                            backgroundColor: "#fff",
                                        }}>
                                            <div dangerouslySetInnerHTML={{ __html: svgString }} />
                                        </Box>
                                    ) : (
                                        <Box>
                                            <Typography>
                                                Nothing to see here.
                                            </Typography>
                                        </Box>
                                    )}
                                </Box>
                            ) : (
                                <Box
                                    // align all to left
                                    sx={{
                                        textAlign: "left"
                                    }}
                                >
                                    <Typography>
                                        {`Category: ${results[selectedResultIndex]["metaData"]["category"]}`}
                                    </Typography>
                                    <Typography>
                                        {`Product: ${results[selectedResultIndex]["metaData"]["product"]}`}
                                    </Typography>
                                    <Typography>
                                        {`Start position: ${results[selectedResultIndex]["metaData"]["startProtocluster"]}`}
                                    </Typography>
                                    <Typography>
                                        {`End position: ${results[selectedResultIndex]["metaData"]["endProtocluster"]}`}
                                    </Typography>
                                    <Typography>
                                        {`Number of motifs: ${results[selectedResultIndex]["query"].length}`}
                                    </Typography>
                                </Box>
                            )}
                        </Box>
                    ) : (
                        <Box 
                            sx={{
                                width: "100%"
                            }}
                        >
                            <Typography 
                                sx={{
                                    textAlign: "left"
                                }}
                            >
                                No result selected.
                            </Typography>
                        </Box>
                    )}
                </Box>
                <Box sx={{ 
                    flexGrow: 1, 
                    overflowY: "auto", 
                    // minHeight: "275px",
                    maxHeight: "275px", 
                    flex: "1 1 50%",
                    backgroundColor: "white",
                    boxShadow: "inset 0 2px 5px rgba(0, 0, 0, 0.1)",
                    border: "1px solid #bbb", 
                    borderRadius: "4px",
                    padding: "10px"
                }}>
                    {results.length ? (
                        <Box>
                            {results.map((result, index) => (
                                <Button
                                    key={index}
                                    variant="contained"
                                    color={selectedResultIndex === index ? "secondary" : "primary"}
                                    onClick={() => handleSelectResultIndex(index)}
                                    sx={{ 
                                        width: "100%", 
                                        margin: "5px 0",
                                        display: "flex",
                                        justifyContent: "flex-start",
                                        alignItems: "center",
                                        textAlign: "left",
                                    }}
                                >
                                    <span style={{ marginRight: "10px" }}>{index + 1}</span>
                                    <span>{result["title"]}</span>
                                </Button>
                            ))}
                        </Box>
                    ) : (
                        <Typography>
                            No results available.
                        </Typography>
                    )}
                </Box>
                
            </Box>
        </Paper>
    );
};

export default ResultsDisplay;