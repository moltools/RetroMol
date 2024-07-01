import React, { useState } from "react";
import { 
    Button, 
    Box, 
    Divider,
    FormControlLabel, 
    IconButton,
    Paper, 
    Radio, 
    RadioGroup, 
    TextField, 
    Typography,
    FormControl,
    InputLabel,
    Select,
    MenuItem 
} from '@mui/material'
import TypographyTooltip from "../common/TypographyTooltip";
import RefreshIcon from '@mui/icons-material/Refresh';
import { toast } from 'react-toastify';

const smilesErythromycin = "CCC1C(C(C(C(=O)C(CC(C(C(C(C(C(=O)O1)C)OC2CC(C(C(O2)C)O)(C)OC)C)OC3C(C(CC(O3)C)N(C)C)O)(C)O)C)C)O)(C)O";
const smilesEnterobactin = "C1C(C(=O)OCC(C(=O)OCC(C(=O)O1)NC(=O)C2=C(C(=CC=C2)O)O)NC(=O)C3=C(C(=CC=C3)O)O)NC(=O)C4=C(C(=CC=C4)O)O";
const smilesDaptomycin = "CCCCCCCCCC(=O)NC(CC1=CNC2=CC=CC=C21)C(=O)NC(CC(=O)N)C(=O)NC(CC(=O)O)C(=O)NC3C(OC(=O)C(NC(=O)C(NC(=O)C(NC(=O)CNC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)CNC3=O)CCCN)CC(=O)O)C)CC(=O)O)CO)C(C)CC(=O)O)CC(=O)C4=CC=CC=C4N)C";
const smilesMicacocidinA = "CCCCCC1=C(C(=CC=C1)[O-])C2=N[C@H](CS2)[C@@H]3N([C@@H](CS3)[C@@H](C(C)(C)C4=N[C@@](CS4)(C)C(=O)[O-])O)C"

const AntismashLogo = ({ style }) => {
    // Make the logo stand out from the background.
    const defaultStyle = { boxShadow: "0 0 10px 5px white" };
  
    return (
        <a href="https://antismash.secondarymetabolites.org/" target="_blank" rel="noreferrer">
            <img 
                src={`${process.env.PUBLIC_URL + "/imgs/antismash_logo.svg"}`} 
                alt="Antismash logo" 
                style={{ ...defaultStyle, ...style }}
            />
        </a>
    );
};

const InputForm = ({ onSubmit }) => {
    const jsonFileNameDefault = "Upload JSON file";

    const [selectedInputType, setSelectedInputType] = useState("smiles");
    const [inputValue, setInputValue] = useState("");
    const [jsonContent, setJsonContent] = useState("");
    const [jsonFileName, setJsonFileName] = useState(jsonFileNameDefault);

    const [selectedExample, setSelectedExample] = useState("");

    const handleChangeExample = (event) => {
        setSelectedExample(event.target.value);
        setSelectedInputType("smiles");
        setInputValue(event.target.value);
    };

    const handleRefresh = () => {
        setInputValue("");
        setJsonContent("");
        setJsonFileName(jsonFileNameDefault);
        setSelectedInputType("smiles");
        setSelectedExample("");
    };

    const handleInputTypeChange = (event) => {
        setInputValue("");
        setJsonContent("");
        setJsonFileName(jsonFileNameDefault);
        setSelectedInputType(event.target.value);
    };

    const handleInputChange = (event) => {
        setInputValue(event.target.value);
    };

    const handleFileChange = (event) => {
        const file = event.target.files[0];
        if (file.size > 20000000) {
            toast.warn("File size exceeds 20 MB. Please upload a smaller file.");
            return;
        }
        setJsonFileName(file.name);
        const reader = new FileReader();
        reader.onload = (event) => {
            setJsonContent(event.target.result);
        };
        reader.readAsText(file);
    };

    const handleFormSubmit = (event) => {
        event.preventDefault();
        if (selectedInputType === "json") {
            onSubmit({
                inputType: selectedInputType,
                inputValue: jsonContent
            });
        } else {
            onSubmit({
                inputType: selectedInputType,
                inputValue: inputValue
            });
        };
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
                minHeight: "385px",
                maxHeight: "385px",
            }}
        >
            <form onSubmit={handleFormSubmit}>
                <Box sx={{ display: "flex", justifyContent: "space-between", alignItems: "center" }}>
                    <Box sx={{ display: "flex", gap: 1 }}>
                        <TypographyTooltip title="AntiSMASH results can be parsed by either passing along a job ID or uploading a result in JSON format.">
                            <Typography variant="h6" gutterBottom>
                                Select input type
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
                <RadioGroup
                    value={selectedInputType}
                    onChange={(event) => handleInputTypeChange(event)}
                    sx={{ mb: 2 }}
                >
                    <FormControlLabel 
                        value="smiles" 
                        control={<Radio />} 
                        label="SMILES" 
                        sx={{ marginLeft: "25px" }}
                    />
                    <Box>
                        <AntismashLogo style={{ width: "25px", transform: "translateY(7px)" }} />
                        <FormControlLabel 
                            value="jobId" 
                            control={<Radio />} 
                            label="Job ID" 
                            sx={{ marginLeft: "0px" }}
                        />
                    </Box>
                    <Box>
                        <AntismashLogo style={{ width: "25px", transform: "translateY(7px)" }} />
                        <FormControlLabel 
                            value="json" 
                            control={<Radio />} 
                            label="JSON" 
                            sx={{ marginLeft: "0px" }}
                        />
                    </Box>
                </RadioGroup>
                {selectedInputType !== "json" && (
                    <TextField
                        sx={{ width: "100%" }}
                        label={selectedInputType === "smiles" ? "Enter SMILES string" : "Enter Job ID"}
                        value={inputValue}
                        onChange={handleInputChange}
                        margin="normal"
                        variant="outlined"
                    />
                )}
                {selectedInputType === "json" &&(
                    <Button
                        variant="outlined"
                        component="label"
                        sx={{ 
                            mt: 3, 
                            width: "100%", 
                            height: "56px", // Align JSON file upload button with text field.
                            transform: "translateY(-8px)" // Align JSON file upload button with text field.
                        }}
                        disabled={jsonContent !== ""}
                    >
                        {jsonFileName}
                        <input
                            type="file"
                            hidden
                            accept=".json"
                            onChange={handleFileChange}
                        />
                    </Button>
                )}
                <Box sx={{ display: "flex", justifyContent: "flex-end", gap: 2 }}>
                    <FormControl sx={{ mt:2, width: 150 }} size="small">
                    <InputLabel id="demo-select-small-label">Example</InputLabel>
                    <Select
                        labelId="demo-select-small-label"
                        id="demo-select-small"
                        value={selectedExample}
                        label="Example"
                        onChange={handleChangeExample}
                    >
                        <MenuItem value="">
                        <em>None</em>
                        </MenuItem>
                        <MenuItem value={smilesDaptomycin}>Daptomycin</MenuItem>
                        {/* <MenuItem value={smilesEnterobactin}>Enterobactin</MenuItem> */}
                        <MenuItem value={smilesErythromycin}>Erythromycin</MenuItem>
                        <MenuItem value={smilesMicacocidinA}>Micacocidin A</MenuItem>
                    </Select>
                    </FormControl>
                    <Button 
                        type="submit" 
                        variant="contained" 
                        sx={{ mt: 2, width: "150px" }}>
                        Submit
                    </Button>
                </Box>
            </form>
        </Paper>
    );
};

export default InputForm;