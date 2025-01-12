import React, { useEffect } from "react";
import { Box } from "@mui/material";
import SmilesDrawer from "smiles-drawer";

class CustomSvgDrawer extends SmilesDrawer.SvgDrawer {
    constructor(options) {
        super(options);
    }

    drawAtomHighlights(highlights) {
        let preprocessor = this.preprocessor;
        let opts = preprocessor.opts;
        let graph = preprocessor.graph;
        let rings = preprocessor.rings;
        let svgWrapper = this.svgWrapper;

        for (var i = 0; i < graph.vertices.length; i++) {
            let vertex = graph.vertices[i];
            let atom = vertex.value;

            for (var j = 0; j < preprocessor.highlight_atoms.length; j++) {
                let highlight = preprocessor.highlight_atoms[j]

                // if atom.bracket !== null, then it is a bracket atom, and we continue
                if (atom.bracket !== null) {
                    if (atom.bracket.isotope === highlight[0]) {
                        svgWrapper.drawAtomHighlight(vertex.position.x, vertex.position.y, highlight[1]);
                    }
                }
            }
        }
    }
}

/**
 * component to draw a molecule from a SMILES string.
 * 
 * @param {number} identifier - Unique identifier for the component.
 * @param {string} smilesStr - SMILES string of the molecule.
 * @returns {React.ReactElement} - The component showing the molecule.
 */ 
const SmilesDrawerContainer = ({ identifier, smilesStr, width, height, highlightAtoms = [] }) => {
    // create a new drawer instance
    let drawer = new CustomSvgDrawer({ width: width, height: height });

    // draw the molecule when the component is mounted
    useEffect(() => {
        let target = `structure-svg-${identifier}`
        let themeName = "light";
        let weights = null;
        let infoOnly = false;
        let weightsNormalized = false;

        SmilesDrawer.parse(smilesStr, function (tree) {
            drawer.draw(tree, target, themeName, weights, infoOnly, highlightAtoms, weightsNormalized);
        });
    }, [smilesStr, highlightAtoms]); // re-draw the molecule when the SMILES string changes

    return (
        <Box key={identifier} sx={{ width: width, height: height }}>
            <svg id={`structure-svg-${identifier}`}/>
        </Box>
    );
};

export default SmilesDrawerContainer;