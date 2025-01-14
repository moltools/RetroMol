import React, { useEffect } from "react";
import { Box } from "@mui/material";
import SmilesDrawer from "smiles-drawer";

class CustomSvgDrawer extends SmilesDrawer.SvgDrawer {
    constructor(options) {
        super(options);
    }
    
    customDrawAtomHighlight(x, y, color = "#03fc9d") {
        let ball = document.createElementNS('http://www.w3.org/2000/svg', 'circle');
        ball.setAttributeNS(null, 'cx', x);
        ball.setAttributeNS(null, 'cy', y);
        ball.setAttributeNS(null, 'r', this.opts.bondLength / 3);
        ball.setAttributeNS(null, 'fill', color);

        this.svgWrapper.highlights.push(ball);
    }

    customDrawBondHighlight(x1, y1, x2, y2, color = "#03fc9d") {
        let line = document.createElementNS('http://www.w3.org/2000/svg', 'line');
        line.setAttributeNS(null, 'x1', x1);
        line.setAttributeNS(null, 'y1', y1);
        line.setAttributeNS(null, 'x2', x2);
        line.setAttributeNS(null, 'y2', y2);
        line.setAttributeNS(null, 'stroke', color);
        line.setAttributeNS(null, 'stroke-width', this.opts.bondLength / 2);

        this.svgWrapper.highlights.push(line);
    }

    drawAtomHighlights(highlights) {
        let preprocessor = this.preprocessor;
        let opts = preprocessor.opts;
        let graph = preprocessor.graph;
        let rings = preprocessor.rings;
        let svgWrapper = this.svgWrapper;

        // highlighted atom ids
        const highlightedAtomIds = [];
        const atomIdToHighlight = {};

        for (var i = 0; i < graph.vertices.length; i++) {
            let vertex = graph.vertices[i];
            let atom = vertex.value;

            for (var j = 0; j < preprocessor.highlight_atoms.length; j++) {
                let highlight = preprocessor.highlight_atoms[j]

                // if atom.bracket !== null, then it is a bracket atom, and we continue
                if (atom.bracket !== null) {
                    if (atom.bracket.isotope === highlight[0]) {
                        this.customDrawAtomHighlight(vertex.position.x, vertex.position.y, highlight[1]);
                        highlightedAtomIds.push(vertex.id);
                        atomIdToHighlight[vertex.id] = highlight[1];
                    }
                }
            }
        }

        // loop over edges
        for (var i = 0; i < graph.edges.length; i++) {
            let edge = graph.edges[i];
            // if edge.sourceId and edge.targetId in highlightedAtomIds, then draw bond highlight, they also need to have same highlight color
            if (highlightedAtomIds.includes(edge.sourceId) && highlightedAtomIds.includes(edge.targetId)) {
                if (atomIdToHighlight[edge.sourceId] === atomIdToHighlight[edge.targetId]) {
                    const sourceVertex = graph.vertices.find(vertex => vertex.id === edge.sourceId); 
                    const targetVertex = graph.vertices.find(vertex => vertex.id === edge.targetId);
                    this.customDrawBondHighlight(sourceVertex.position.x, sourceVertex.position.y, targetVertex.position.x, targetVertex.position.y, atomIdToHighlight[edge.sourceId]);
                }
            }
        }

        // loop over all atoms and set atom.bracket to null
        for (var i = 0; i < graph.vertices.length; i++) {
            let vertex = graph.vertices[i];
            let atom = vertex.value;
            atom.bracket = null;

            // make sure COOH is drawn fully instead of displayed with text
            if (atom.element === "C") {
                if (atom.hasAttachedPseudoElements) {
                }
                atom.hasAttachedPseudoElements = false;
            }

            if (atom.element === "O") {
                atom.isDrawn = true;
            }
        }

        // loop over all bonds
        for (var i = 0; i < graph.edges.length; i++) {
            let edge = graph.edges[i];
            // if aromatic bond
            if (edge.isPartOfAromaticRing) {
                // set to false to prevent drawing of double bond rings in ring
                edge.isPartOfAromaticRing = false;
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