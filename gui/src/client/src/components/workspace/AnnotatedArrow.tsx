import React from "react";
import Box from "@mui/material/Box";
import Typography from "@mui/material/Typography";

// A labeled "->" used between stages of a structure/readout diagram (e.g.
// "structure -> Linearization -> primary sequence" for a compound, or
// "genes/domains -> Readout -> primary sequence" for a gene cluster).
export function AnnotatedArrow({ annotation }: { annotation: string }) {
  return (
    <Box
      sx={{
        display: 'flex',
        flexDirection: 'column',
        alignItems: 'center',
        justifyContent: 'center',
        flexShrink: 0, // don’t let it compress
      }}
    >
      <Typography variant='caption' gutterBottom>
        {annotation}
      </Typography>
      <Box
        component='svg'
        width={140} // total width of arrow+shaft
        height={32}
        viewBox='0 0 140 32'
      >
        {/* horizontal line */}
        <line
          x1='0'
          y1='16'
          x2='120'
          y2='16'
          stroke='currentColor'
          strokeWidth='2'
        />
        {/* arrow head */}
        <polygon
          points='120,8 120,24 140,16'
          fill='currentColor'
        />
      </Box>
    </Box>
  );
}
