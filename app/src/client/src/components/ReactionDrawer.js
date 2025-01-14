import React, { useEffect } from "react";
import { Box } from "@mui/material";
import SmilesDrawer from "smiles-drawer";

const ReactionDrawerContainer = ({ identifier, smilesStr, width, height }) => {
    // create a new drawer instance
    // let drawer = new SmilesDrawer.SvgDrawer({ width: width, height: height });
    let drawer = new SmilesDrawer.ReactionDrawer(
        { width: width, height: height }, 
        { scale: 1.0 }
    );

    // draw the molecule when the component is mounted
    useEffect(() => {
        let target = `reaction-svg-${identifier}`
        let themeName = "light";

        SmilesDrawer.parseReaction(smilesStr, function (tree) {
            console.log(drawer);
            drawer.draw(tree, target, themeName);
        });
    }, [smilesStr]); // re-draw the molecule when the SMILES string changes

    return (
        <Box key={identifier} sx={{ width: width, height: height }}>
            <svg id={`reaction-svg-${identifier}`}/>
        </Box>
    );
};

export default ReactionDrawerContainer;