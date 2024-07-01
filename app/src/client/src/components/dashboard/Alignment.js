import React, { useState, useRef } from "react";
import { Slider, Typography, Box,
    Button, Table, TableBody, TableCell, TableContainer, TableHead, TableRow, Paper, Tooltip } from "@mui/material";
import * as htmlToImage from 'html-to-image';
import { toast } from "react-toastify";

import { ZoomIn, ZoomOut } from "@mui/icons-material";


const Alignment = ({ data }) => {
    const [zoom, setZoom] = useState(1);
    const tableRef = useRef(null);
    const sequenceLength = data[0].motifCode.length;

    const handleZoomChange = (event, newValue) => {
        setZoom(newValue);
    };

    // Define colors based on starting characters.
    const getColor = (motif) => {
        const motifType = motif.motifType;
        const polyketideType = motif.polyketideType;

        if (motifType === "peptide") {
            return "#ed9ea8";
        } else if (motifType === "polyketide" && polyketideType === "A") {
            return "#209bef";
        } else if (motifType === "polyketide" && polyketideType === "B") {
            return "#fede57";
        } else if (motifType === "polyketide" && polyketideType === "C") {
            return "#47c774";
        } else if (motifType === "polyketide" && polyketideType === "D") {
            return "#ff3960";
        } else if (motifType === "gap") {
            return "transparent";
        } else {
            return "#ccc";
        }
    };

    // Define labels based on starting characters.
    const getLabel = (motif) => {
        const motifType = motif.motifType;
        var polyketideType = motif.polyketideType;
        var polyketideDecor = motif.polyketideDecor;
        var cid = motif.peptideCid;
        
        // Set Any to wildcard.
        if (motifType === "polyketide" && polyketideType === "Any") {
            polyketideType = "*";
        };
        if (motifType === "polyketide" && polyketideDecor === "Any") {
            polyketideDecor = "*";
        }
        if (motifType === "peptide" && cid === "Any") {
            cid = "*";
        }

        // Decide on label.
        if (motifType === "peptide") {
            return cid;
        } else if (motifType === "polyketide") {
            return polyketideType + polyketideDecor;
        } else if (motifType === "gap") {
            return "";
        } else {
            return "?";
        }
    };

    const handleDownloadPng = () => {
        htmlToImage.toPng(tableRef.current)
            .then(function (dataUrl) {
                var link = document.createElement('a');
                link.download = 'alignment.png';
                link.href = dataUrl;
                link.click();
            })
            .catch(function (error) {
                toast.error('Error generating image:', error);
            });
    };

    return (
        <div>
            <Box display="flex" flexDirection="column" alignItems="left" marginBottom="1rem">
                <Box display="flex" alignItems="center" gap={1}>
                    <ZoomOut />
                    <Slider
                        value={zoom}
                        min={0.5}
                        max={2}
                        step={0.1}
                        onChange={handleZoomChange}
                        aria-labelledby="zoom-slider"
                        style={{ width: '200px' }}
                    />
                    <ZoomIn />
                </Box>
                {/* <Button
                    variant="contained"
                    onClick={handleDownloadPng}
                    style={{ marginBottom: "1rem", marginLeft: "1rem" }}
                >
                    Download as PNG
                </Button> */}
            </Box>
            <div ref={tableRef} style={{ overflowX: "auto", overflowY: "auto", transform: `scale(${zoom})`, transformOrigin: '0 0' }}>
                <TableContainer component={Paper}>
                    <Table style={{ borderCollapse: "collapse", width: "auto" }}>
                        <TableHead>
                            <TableRow>
                                <TableCell>Identifier</TableCell>
                                {Array.from(Array(sequenceLength).keys()).map(index => (
                                    <TableCell key={index} style={{ textAlign: "center", padding: 0}}>
                                        {index + 1}
                                    </TableCell>
                                ))}
                                <TableCell style={{ textAlign: "left" }}>Bioactivity</TableCell>
                                <TableCell style={{ textAlign: "left" }}>Genus</TableCell>
                            </TableRow>
                        </TableHead>
                        <TableBody>
                            {data.map((item, index) => (
                                <TableRow key={index}>
                                    <TableCell 
                                        style={{ 
                                            width: "150px",
                                            minWidth: "150px",
                                            maxWidth: "150px",
                                            height: "50px",
                                            minHeight: "50px",
                                            maxHeight: "50px",
                                            overflow: "hidden",
                                            textOverflow: "ellipsis",
                                            whiteSpace: "nowrap",
                                        }}
                                    >
                                        <Tooltip title={item.motifCodeIdentifier} placement="top" arrow>
                                            {item.motifCodeIdentifier.startsWith("NPA") ? (
                                                <a
                                                    href={`https://www.npatlas.org/explore/compounds/${item.motifCodeIdentifier}`}
                                                    target="_blank"
                                                    rel="noreferrer"
                                                >
                                                    {item.motifCodeIdentifier}
                                                </a>
                                            ) : (
                                                item.motifCodeIdentifier
                                            )}
                                        </Tooltip>
                                    </TableCell>
                                    {item.motifCode.map((motif, index) => (
                                        <TableCell
                                            key={index}
                                            style={{
                                                backgroundColor: getColor(motif),
                                                border: "1px solid #ccc",
                                                width: "50px",
                                                minWidth: "50px",
                                                maxWidth: "50px",
                                                height: "50px",
                                                minHeight: "50px",
                                                maxHeight: "50px",
                                                overflow: "hidden",
                                                textOverflow: "ellipsis",
                                                whiteSpace: "nowrap",
                                                paddingRight: "5px",
                                                paddingLeft: "5px",
                                                textAlign: "center"
                                            }}
                                        >   
                                            <Tooltip title={getLabel(motif)} placement="top" arrow>
                                                <span>{getLabel(motif)}</span>
                                            </Tooltip>
                                        </TableCell>
                                    ))}
                                    <Tooltip title={item.bioactivities.length ? item.bioactivities.join(", ") : "N/A"} placement="top" arrow>
                                        <TableCell style={{ 
                                            textAlign: "left", 
                                            whiteSpace: "nowrap",
                                            overflow: "hidden", 
                                            textOverflow: "ellipsis", 
                                            maxWidth: 300 // Adjust the width as needed 
                                        }}>
                                            {item.bioactivities.length ? item.bioactivities.join(", ").charAt(0).toUpperCase() + item.bioactivities.join(", ").slice(1) : "N/A"}
                                        </TableCell>
                                    </Tooltip>
                                    <Tooltip title={item.genus.length ? item.genus.join(", ") : "N/A"} placement="top" arrow>
                                        <TableCell style={{ 
                                            textAlign: "left", 
                                            whiteSpace: "nowrap",
                                            overflow: "hidden",
                                            textOverflow: "ellipsis",
                                            maxWidth: 300 // Adjust the width as needed 
                                        }}>
                                            {item.genus.length ? item.genus.join(", ").charAt(0).toUpperCase() + item.genus.join(", ").slice(1) : "N/A"}
                                        </TableCell>
                                    </Tooltip>
                                </TableRow>
                            ))}
                        </TableBody>
                    </Table>
                </TableContainer>
            </div>
        </div>
    );
};

export default Alignment;