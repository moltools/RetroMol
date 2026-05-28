import React from "react";
import Box from "@mui/material/Box";
import Stack from "@mui/material/Stack";
import Typography from "@mui/material/Typography";
import Tooltip from "@mui/material/Tooltip";
import { CSS } from "@dnd-kit/utilities";
import { useSortable } from "@dnd-kit/sortable";
import DragIndicatorIcon from "@mui/icons-material/DragIndicator";
import InfoOutlineIcon from "@mui/icons-material/InfoOutline";
import { MsaRow } from "./query_types";
import { DialogRowInfo } from "./DialogRowInfo";

interface SortableRowProps {
  row: MsaRow;
  labelWidth: number;
  columnTemplate: string;
  children: React.ReactNode;
  hasRowInfo: boolean;
};

const fmt = (v: number | null, digits: number) =>
  v == null || Number.isNaN(v) ? "" : v.toFixed(digits);

export const SortableRow: React.FC<SortableRowProps> = ({
  row,
  labelWidth,
  columnTemplate,
  children,
  hasRowInfo,
}) => {
  const {
    attributes,
    listeners,
    setNodeRef,
    setActivatorNodeRef,
    transform,
    transition,
  } = useSortable({ id: row.id });

  const alignText = fmt(row.alignment_score, 2);
  const cosineText = fmt(row.cosine_score, 2);
  const matchText = fmt(row.match_score, 2);
  const scoreBlockWidth = 40;

  const [openRowInfo, setOpenRowInfo] = React.useState(false);

  const handleRowInfo = (event: React.MouseEvent<SVGSVGElement>) => {
    setOpenRowInfo(true);
  };

  return (
    <>
      <Box
        ref={setNodeRef}
        sx={{
          display: "grid",
          gridTemplateColumns: columnTemplate,
          columnGap: 1,
          alignItems: "center",
          transform: transform ? CSS.Transform.toString(transform) : undefined,
          transition,
        }}
      >
        <Box
          sx={{
            m: 0,
            p: 0,
            height: 20,
            fontWeight: 600,
            zIndex: 100,
            display: "flex",
            alignItems: "center",
            justifyContent: "space-between",
            px: 1,
            borderRadius: 1,
            backgroundColor: "background.paper",
            border: "1px solid transparent",
            borderColor: "transparent",
          }}
        >
          <Box
            ref={setActivatorNodeRef}
            {...listeners}
            {...attributes}
            sx={{
              display: "flex",
              alignItems: "center",
              cursor: "grab",
            }}
            onClick={(e) => e.stopPropagation()} // don't trigger center selection on drag
          >
            <DragIndicatorIcon fontSize="small" />
          </Box>
          {hasRowInfo && (
            <Box>
              <InfoOutlineIcon
                fontSize="small"
                onClick={handleRowInfo}
                sx={{
                  mt: "5px",
                  cursor: "pointer",
                  color: "text.secondary",
                  "&:hover": { color: "text.primary" },
                }}
              />
            </Box>
          )}
          <Stack
            direction="row"
            alignItems="center"
            justifyContent="space-between"
            sx={{ width: "100%", pl: 1 }}
          >
            <Tooltip title={row.name || row.id} arrow>
              <Typography
                variant="body2"
                component="span"
                sx={{
                  fontWeight: 600,
                  maxWidth: labelWidth - 100,
                  lineHeight: "20px",
                  overflow: "hidden",
                  textOverflow: "ellipsis",
                  whiteSpace: "nowrap",
                  zIndex: 101,
                  userSelect: "none",
                }}
              >
                {row.name || row.id}
              </Typography>
            </Tooltip>

            {/* Scores */}
            <Stack
              direction="column"
              sx={{
                width: scoreBlockWidth,
                flex: `0 0 ${scoreBlockWidth}px`,
                alignItems: "flex-end",
                justifyContent: "center",
                lineHeight: 1,
              }}
            >
              <Tooltip title={`Alignment score: ${alignText}`} placement="right" arrow>
                <Typography
                  variant="caption"
                  sx={{
                    lineHeight: 1,
                    fontFamily: "monospace",
                    fontSize: "0.70rem",
                  }}
                >
                  {alignText}
                </Typography>
              </Tooltip>
              <Tooltip title={`Cosine score: ${cosineText}`} placement="right" arrow>
                <Typography
                  variant="caption"
                  sx={{
                    lineHeight: 1,
                    fontFamily: "monospace",
                    fontSize: "0.70rem",
                  }}
                >
                  {cosineText}
                </Typography>
              </Tooltip>
              <Tooltip title={`Ratio of items aligned: ${matchText}`} placement="right" arrow>
                <Typography
                  variant="caption"
                  sx={{
                    lineHeight: 1,
                    fontFamily: "monospace",
                    fontSize: "0.70rem",
                  }}
                >
                  {matchText}
                </Typography>
              </Tooltip>
            </Stack>
          </Stack>
        </Box>

        {/* Motifs */}
        {children}
      </Box>

      <DialogRowInfo
        open={openRowInfo}
        onClose={() => setOpenRowInfo(false)}
        msaRow={row}
      />
    </>
  );
};