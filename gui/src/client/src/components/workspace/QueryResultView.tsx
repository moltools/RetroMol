import React from "react";
import Typography from "@mui/material/Typography";
import Box from "@mui/material/Box";
import Chip from "@mui/material/Chip";
import Stack from "@mui/material/Stack";
import Tooltip from "@mui/material/Tooltip";

import { DndContext, DragEndEvent } from "@dnd-kit/core";
import { SortableContext, arrayMove, verticalListSortingStrategy } from "@dnd-kit/sortable";

import { SortableRow } from "./SortableRow";
import { useNotifications } from "../NotificationProvider";

import type { QueryResult, MsaRow } from "./query_types";
import { PROTECTED_NAME_TO_CODE, isPolyketideMotif, makeToDisplayName, renderChipLabel, renderTooltipLabel } from "./query_motif";
import { getMotifColor } from "./query_colors";
import { buildMsaSvg, sequenceLength } from "./query_svg";
import { QueryResultToolbar } from "./QueryResultToolbar";

type QueryResultViewProps = {
  result: QueryResult;
};

export const QueryResultView: React.FC<QueryResultViewProps> = ({ result }) => {
  const { pushNotification } = useNotifications();

  const [msa, setMsa] = React.useState<MsaRow[]>(result.msa);

  React.useEffect(() => {
    setMsa(result.msa);
  }, [result]);

  // stable name->code mapping
  const toDisplayName = React.useMemo(
    () => makeToDisplayName(PROTECTED_NAME_TO_CODE),
    []
  );

  const msaLength = React.useMemo(() => {
    if (msa.length === 0) return 0;
    return Math.max(...msa.map((r) => sequenceLength(r.sequence)));
  }, [msa]);

  // Zoom
  const [zoom, setZoom] = React.useState<number>(1.0);
  const handleZoomIn = () => setZoom((z) => Math.min(z + 0.1, 3.0));
  const handleZoomOut = () => setZoom((z) => Math.max(z - 0.1, 0.5));
  const handleZoomReset = () => setZoom(1.0);

  const invertMsaMotifOrder = React.useCallback(() => {
    setMsa((prev) =>
      prev.map((row) => ({
        ...row,
        sequence: [...row.sequence]
          .reverse()
          .map((seq) => ({ ...seq, sequence: [...seq.sequence].reverse() })),
      }))
    );
  }, []);

  const labelWidth = 250;
  const motifWidth = 50 * zoom;
  const colTemplate = React.useMemo(() => {
    return `${labelWidth}px repeat(${msaLength}, ${motifWidth}px) 1fr`;
  }, [labelWidth, msaLength, motifWidth]);

  const handleDragEnd = React.useCallback((event: DragEndEvent) => {
    const { active, over } = event;
    if (!over) return;

    const activeId = active.id as string;
    const overId = over.id as string;
    if (activeId === overId) return;

    setMsa((prev) => {
      const rowIds = new Set(prev.map((r) => r.id));
      const isRowDrag = rowIds.has(activeId) && rowIds.has(overId);

      if (isRowDrag) {
        const fromIndex = prev.findIndex((r) => r.id === activeId);
        const toIndex = prev.findIndex((r) => r.id === overId);
        if (fromIndex === -1 || toIndex === -1) return prev;
        return arrayMove(prev, fromIndex, toIndex);
      }

      let rowIndex = -1;
      let fromCol = -1;
      let toCol = -1;

      for (let r = 0; r < prev.length; r++) {
        const seq = prev[r].sequence;
        const aIdx = seq.findIndex((s) => s.id === activeId);
        const oIdx = seq.findIndex((s) => s.id === overId);
        if (aIdx !== -1 && oIdx !== -1) {
          rowIndex = r;
          fromCol = aIdx;
          toCol = oIdx;
          break;
        }
      }

      if (rowIndex === -1) return prev;

      const row = prev[rowIndex];
      const newSeq = arrayMove([...row.sequence], fromCol, toCol);
      return prev.map((r, idx) => (idx === rowIndex ? { ...r, sequence: newSeq } : r));
    });
  }, []);

  const handleDownloadMsaSvg = React.useCallback(() => {
    const svg = buildMsaSvg({ msa, msaLength, toDisplayName });

    if (!svg) {
      pushNotification("No visible sequences to download.", "warning");
      return;
    }

    const blob = new Blob([svg], { type: "image/svg+xml;charset=utf-8" });
    const url = URL.createObjectURL(blob);
    const link = document.createElement("a");
    link.href = url;
    link.download = "msa.svg";
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
    URL.revokeObjectURL(url);
  }, [msa, msaLength, pushNotification, toDisplayName]);

  return (
    <div>
      <Typography component="h1" variant="subtitle1">
        Query results
      </Typography>

      <QueryResultToolbar
        onZoomOut={handleZoomOut}
        onZoomIn={handleZoomIn}
        onZoomReset={handleZoomReset}
        onInvertMotifs={invertMsaMotifOrder}
        onDownloadSvg={handleDownloadMsaSvg}
      />

      <Stack direction="column" spacing={1}>
        <Box sx={{ width: "100%", overflowX: "auto", overflowY: "hidden", py: 2 }}>
          <DndContext onDragEnd={handleDragEnd}>
            <SortableContext items={msa.map((item) => item.id)} strategy={verticalListSortingStrategy}>
              <Box sx={{ display: "flex", flexDirection: "column", width: "100%", gap: 3.5, py: 1 }}>
                {msa.map((row, rowIndex) => (
                  <Box key={row.id} sx={{ position: "relative", width: "fit-content" }}>
                    <SortableRow
                      key={row.id}
                      row={row}
                      labelWidth={labelWidth}
                      columnTemplate={colTemplate}
                      hasRowInfo={rowIndex > 0}
                    >
                      {row.sequence.map((subseq) => {
                        const allGaps = subseq.sequence.every((it) => it.isGap);

                        return (
                          <Box
                            key={subseq.id}
                            sx={{
                              display: "flex",
                              flexDirection: "row",
                              gridColumn: `span ${subseq.sequence.length}`,
                              position: "relative",
                              zIndex: 110,
                              visibility: allGaps ? "hidden" : "visible",
                              pointerEvents: allGaps ? "none" : "auto",
                              minHeight: 24,
                              "&::before": {
                                content: '""',
                                position: "absolute",
                                left: 0,
                                right: 0,
                                top: -20,
                                height: "32px",
                                borderRadius: "4px",
                                backgroundColor: "divider",
                                pointerEvents: "none",
                                opacity: allGaps ? 0 : 1,
                              },
                            }}
                          >
                            <Typography
                              variant="caption"
                              sx={{
                                position: "absolute",
                                top: -20,
                                left: 4,
                                maxWidth: `${subseq.sequence.length * motifWidth - 8}px`,
                                whiteSpace: "nowrap",
                                overflow: "hidden",
                                textOverflow: "ellipsis",
                                fontFamily: "monospace",
                                fontSize: "0.65rem",
                                color: "text.primary",
                                pointerEvents: "none",
                              }}
                            >
                              {subseq.name || subseq.id}
                            </Typography>

                            <Box sx={{ display: "flex", flexDirection: "row", gap: 1 }}>
                              {subseq.sequence.map((item) =>
                                item.isGap ? (
                                  <Box
                                    key={item.id}
                                    sx={{
                                      m: 0,
                                      width: motifWidth,
                                      height: 20,
                                      backgroundColor: "transparent",
                                      display: "flex",
                                      alignItems: "center",
                                      justifyContent: "center",
                                      mt: "2px",
                                      zIndex: 1000,
                                    }}
                                  >
                                    <Box
                                      sx={{
                                        width: Math.max(5, 3 * zoom),
                                        height: Math.max(5, 3 * zoom),
                                        borderRadius: "50%",
                                        backgroundColor: "text.secondary",
                                      }}
                                    />
                                  </Box>
                                ) : (
                                  <Box
                                    key={item.id}
                                    component="span"
                                    sx={{
                                      backgroundColor: "background.paper",
                                      borderRadius: "10px",
                                      display: "inline-block",
                                      zIndex: 99,
                                    }}
                                  >
                                    <Tooltip
                                      title={
                                        <span>
                                          Motif name is{" "}
                                          {isPolyketideMotif(item.name)
                                            ? renderChipLabel(item.name, toDisplayName)
                                            : renderTooltipLabel(item.name, toDisplayName)}
                                        </span>
                                      }
                                      arrow
                                    >
                                      <Chip
                                        label={renderChipLabel(item.name, toDisplayName)}
                                        size="small"
                                        sx={{
                                          borderRadius: "10px",
                                          width: motifWidth,
                                          height: 20,
                                          textAlign: "center",
                                          backgroundColor: getMotifColor(item.name || "") || "background.paper",
                                          border: "1px solid",
                                          borderColor: getMotifColor(item.name || "") || "primary.main",
                                          zIndex: 100,
                                          fontSize: `${Math.max(0.5, Math.min(1.0, zoom)) * 0.75}rem`,
                                          transition: "background-color 0s ease-in-out",
                                          mt: "-2px",
                                        }}
                                      />
                                    </Tooltip>
                                  </Box>
                                )
                              )}
                            </Box>
                          </Box>
                        );
                      })}
                    </SortableRow>
                  </Box>
                ))}
              </Box>
            </SortableContext>
          </DndContext>
        </Box>
      </Stack>
    </div>
  );
};