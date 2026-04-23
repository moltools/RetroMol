import React, { useEffect } from 'react';
import SmilesDrawer from 'smiles-drawer';
import { Box } from '@mui/material';
import { useColorScheme } from '@mui/material/styles';

class CustomSvgDrawer extends SmilesDrawer.SvgDrawer {
    constructor(options) {
        super(options);

        const themeOverrides = {
            light: {
                C: '#000000',
                O: '#000000',
                N: '#000000',
                S: '#000000',
                H: '#000000',
            },
            dark: {
                C: '#ffffff',
                O: '#ffffff',
                N: '#ffffff',
                S: '#ffffff',
                H: '#ffffff',
            }
        };

        Object.entries(themeOverrides).forEach(([themeName, overrides]) => {
            this.opts.themes[themeName] = {
                // copy whatever was in the built-in theme
                ...(this.opts.themes[themeName] || {}),
                // then overwrite with your custom colors
                ...overrides
            };
        });
    };

    customDrawAtomHighlight(x, y, color = '#03fc9d') {
        let ball = document.createElementNS('http://www.w3.org/2000/svg', 'circle');
        ball.setAttributeNS(null, 'cx', x);
        ball.setAttributeNS(null, 'cy', y);
        ball.setAttributeNS(null, 'r', this.opts.bondLength / 3);
        ball.setAttributeNS(null, 'fill', color);

        this.svgWrapper.highlights.push(ball);
    };

    customDrawBondHighlight(x1, y1, x2, y2, color = '#03fc9d') {
        let line = document.createElementNS('http://www.w3.org/2000/svg', 'line');
        line.setAttributeNS(null, 'x1', x1);
        line.setAttributeNS(null, 'y1', y1);
        line.setAttributeNS(null, 'x2', x2);
        line.setAttributeNS(null, 'y2', y2);
        line.setAttributeNS(null, 'stroke', color);
        line.setAttributeNS(null, 'stroke-width', this.opts.bondLength / 2);

        this.svgWrapper.highlights.push(line);
    };

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
                    };
                };
            };
        };

        // loop over edges
        for (var i = 0; i < graph.edges.length; i++) {
            let edge = graph.edges[i];
            // if edge.sourceId and edge.targetId in highlightedAtomIds, then draw bond highlight, they also need to have same highlight color
            if (highlightedAtomIds.includes(edge.sourceId) && highlightedAtomIds.includes(edge.targetId)) {
                if (atomIdToHighlight[edge.sourceId] === atomIdToHighlight[edge.targetId]) {
                    const sourceVertex = graph.vertices.find(vertex => vertex.id === edge.sourceId);
                    const targetVertex = graph.vertices.find(vertex => vertex.id === edge.targetId);
                    this.customDrawBondHighlight(sourceVertex.position.x, sourceVertex.position.y, targetVertex.position.x, targetVertex.position.y, atomIdToHighlight[edge.sourceId]);
                };
            };
        };

        // loop over all atoms and set atom.bracket to null
        for (var i = 0; i < graph.vertices.length; i++) {
            let vertex = graph.vertices[i];
            let atom = vertex.value;
            atom.bracket = null;

            // make sure COOH is drawn fully instead of displayed with text
            if (atom.element === 'C') {
                if (atom.hasAttachedPseudoElements) {
                }
                atom.hasAttachedPseudoElements = false;
            };

            if (atom.element === 'O') {
                atom.isDrawn = true;
            };

            if (atom.element === '*') {  // always draw wildcards
                atom.isDrawn = true;
            };
        };

        // loop over all bonds
        for (var i = 0; i < graph.edges.length; i++) {
            let edge = graph.edges[i];
            // if aromatic bond
            if (edge.isPartOfAromaticRing) {
                // set to false to prevent drawing of double bond rings in ring
                edge.isPartOfAromaticRing = false;
            };
        };

    };
};

/**
 * @typedef {[number, string]} HighlightAtom
 * @typedef {Object} Props
 * @property {string} identifier                    unique ID for the SVG element
 * @property {string} smiles                        SMILES string to draw
 * @property {number} size                          width & height of the drawing
 * @property {HighlightAtom[]} [highlightAtoms]     array of `[atomNumber, color]`
 * @property {string} [themeOverride]               force “light” or “dark” drawing theme
 */

/**
 * @param {Props} props
 */
const SmilesDrawerContainer = ({ identifier, smiles, size, highlightAtoms = [], themeOverride = '' }) => {
    // create a new drawer instance
    let drawer = new CustomSvgDrawer({ width: size, height: size });

    const { mode, systemMode, setMode } = useColorScheme();

    // draw the molecule when the component is mounted
    useEffect(() => {
        let target = `structure-svg-${identifier}`
        let themeName = themeOverride !== '' ? themeOverride : (systemMode !== undefined) ? systemMode : mode;
        let weights = null;
        let infoOnly = false;
        let weightsNormalized = false;

        SmilesDrawer.parse(smiles, function (tree) {
            drawer.draw(tree, target, themeName, weights, infoOnly, highlightAtoms, weightsNormalized);
        });
    }, [smiles, highlightAtoms, size]);

    return (
        <Box key={identifier} sx={{ width: size, height: size }}>
            <svg id={`structure-svg-${identifier}`}/>
        </Box>
    );
};

export default SmilesDrawerContainer;