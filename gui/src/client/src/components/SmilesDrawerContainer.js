import React, { useEffect, useState } from 'react';
import SmilesDrawer from 'smiles-drawer';
import { Box, Typography } from '@mui/material';
import { useColorScheme } from '@mui/material/styles';


class CustomSvgDrawer extends SmilesDrawer.SvgDrawer {
    constructor(options, showIsotopes = false, orientationTags = null) {
        super(options);
        this.showIsotopes = showIsotopes;
        // { startTags, endTags }: isotope tags (see retromol.chem.tagging) of the
        // first and last blocks in a generated sequence -- when given, the drawing
        // is mirrored horizontally (if needed) so the start block ends up on the
        // left, matching reading order of the primary sequence shown alongside it.
        this.orientationTags = orientationTags;

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

    prepareTinAsRightWildcard(element = 'Sn') {
        const graph = this.preprocessor?.graph;
        if (!graph) return;

        const vertices = graph.vertices.filter((v) => v.position);
        if (!vertices.length) return;

        const snVertex = vertices.find((v) => v.value?.element === element);
        if (!snVertex) return;

        // 1) orient so Sn is on the right
        const xs = vertices.map((v) => v.position.x);
        const centerX = (Math.min(...xs) + Math.max(...xs)) / 2;

        if (snVertex.position.x < centerX) {
            for (const vertex of vertices) {
                vertex.position.x = 2 * centerX - vertex.position.x;
            }
            for (const ring of this.preprocessor?.rings ?? []) {
                ring.center.x = 2 * centerX - ring.center.x;
            }
        }

        // 2) replace Sn with wildcard atom
        const atom = snVertex.value;
        atom.element = 'H';
        atom.isDrawn = false;
        atom.hasAttachedPseudoElements = false;
        atom.bracket = null;
    }

    orientSequenceLeftToRight(startTags, endTags) {
        const graph = this.preprocessor?.graph;
        if (!graph || !startTags?.length || !endTags?.length) return;

        const vertices = graph.vertices.filter((v) => v.position);
        if (!vertices.length) return;

        const avgX = (tags) => {
            const matched = vertices.filter((v) => v.value?.bracket && tags.includes(v.value.bracket.isotope));
            if (!matched.length) return null;
            return matched.reduce((sum, v) => sum + v.position.x, 0) / matched.length;
        };

        const startX = avgX(startTags);
        const endX = avgX(endTags);
        if (startX === null || endX === null || startX <= endX) return;

        const xs = vertices.map((v) => v.position.x);
        const centerX = (Math.min(...xs) + Math.max(...xs)) / 2;
        for (const vertex of vertices) {
            vertex.position.x = 2 * centerX - vertex.position.x;
        }
        for (const ring of this.preprocessor?.rings ?? []) {
            ring.center.x = 2 * centerX - ring.center.x;
        }
    }

    drawAtomHighlights(highlights) {
        this.prepareTinAsRightWildcard('Sn');
        if (this.orientationTags) {
            this.orientSequenceLeftToRight(this.orientationTags.startTags, this.orientationTags.endTags);
        }

        let preprocessor = this.preprocessor;
        let graph = preprocessor.graph;

        // highlighted atom ids
        const highlightedAtomIds = [];
        const atomIdToHighlight = {};

        for (let i = 0; i < graph.vertices.length; i++) {
            let vertex = graph.vertices[i];
            let atom = vertex.value;

            for (let j = 0; j < preprocessor.highlight_atoms.length; j++) {
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
        for (let i = 0; i < graph.edges.length; i++) {
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

        // loop over all atoms and, unless isotopes were asked to be kept (used to show
        // R-group numbers, e.g. [1*]/[2*], in the rules browser and motif hover card),
        // clear atom.bracket so isotope/charge annotations don't clutter the drawing.
        for (let i = 0; i < graph.vertices.length; i++) {
            let vertex = graph.vertices[i];
            let atom = vertex.value;
            if (!this.showIsotopes) {
                atom.bracket = null;
            }

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
        for (let i = 0; i < graph.edges.length; i++) {
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
 * @property {boolean} [showIsotopes]               draw isotope numbers (e.g. R-group
 *                                                   tags like [1*]/[2*]) instead of hiding them
 * @property {{startTags: number[], endTags: number[]}} [orientationTags]
 *                                                   isotope tags of the first/last blocks in a
 *                                                   generated sequence -- mirrors the drawing
 *                                                   horizontally (if needed) so the start block
 *                                                   ends up on the left
 */

/**
 * @param {Props} props
 */
const SmilesDrawerContainer = ({ identifier, smiles, size, highlightAtoms = [], themeOverride = '', showIsotopes = false, orientationTags = null }) => {
    const { mode, systemMode } = useColorScheme();
    const [error, setError] = useState(null);

    // draw the molecule when the component is mounted
    useEffect(() => {
        setError(null);

        let target = `structure-svg-${identifier}`
        let themeName = themeOverride !== '' ? themeOverride : (systemMode !== undefined) ? systemMode : mode;
        let weights = null;
        let infoOnly = false;
        let weightsNormalized = false;

        // A fresh drawer instance per draw call, since it isn't safe to reuse
        // once it holds a graph for a previous (possibly differently-sized) SMILES.
        //
        // padding is bumped above the library default (10) because that default
        // only leaves room for atom labels -- it doesn't account for the extra
        // radius a highlight circle (customDrawAtomHighlight, r = bondLength / 3)
        // draws around a highlighted atom, so a highlighted atom sitting at the
        // very edge of the layout (as the first/last block often does once
        // orientationTags pins it there) gets its highlight clipped by the SVG's
        // own viewBox.
        let drawer = new CustomSvgDrawer({ width: size, height: size, padding: 24 }, showIsotopes, orientationTags);

        try {
            SmilesDrawer.parse(
                smiles,
                function (tree) {
                    try {
                        drawer.draw(tree, target, themeName, weights, infoOnly, highlightAtoms, weightsNormalized);
                    } catch (drawErr) {
                        console.error('SmilesDrawerContainer: draw failed', drawErr);
                        setError('Could not render this structure.');
                    }
                },
                function (parseErr) {
                    console.error('SmilesDrawerContainer: parse failed', parseErr);
                    setError('Could not parse this SMILES.');
                }
            );
        } catch (err) {
            // Some malformed input can throw synchronously instead of hitting the error callback
            console.error('SmilesDrawerContainer: unexpected error', err);
            setError('Could not render this structure.');
        }
    }, [identifier, smiles, highlightAtoms, size, themeOverride, showIsotopes, orientationTags, mode, systemMode]);

    if (error) {
        return (
            <Box key={identifier} sx={{ width: size, height: size, display: 'flex', alignItems: 'center', justifyContent: 'center', textAlign: 'center', p: 1 }}>
                <Typography variant="caption" color="text.secondary">{error}</Typography>
            </Box>
        );
    }

    return (
        <Box key={identifier} sx={{ width: size, height: size }}>
            <svg id={`structure-svg-${identifier}`}/>
        </Box>
    );
};

export default SmilesDrawerContainer;