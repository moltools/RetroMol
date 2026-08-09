// Renders a saved discovery query's results -- either a live session item (queued/
// processing/done/error, updated over SSE like DialogViewItem) or an uploaded result
// file loaded purely client-side (see WorkspaceDiscovery's "Load result file"). Both
// cases pass the same DiscoveryQueryItem shape, so most of this component doesn't
// need to know which one it's looking at -- the one exception is `payload`: getSession
// strips it from every item (same as it already does for compounds/clusters), so a
// live item's payload has to be fetched separately once it's done, via
// getDiscoveryQueryResult. An uploaded item was serialized *before* that stripping
// (straight from this same component's own resolved payload, see handleDownload in
// WorkspaceDiscovery), so it already carries its payload and needs no fetch --
// pass sessionId={null} for that case.

import React from "react";
import Alert from "@mui/material/Alert";
import Box from "@mui/material/Box";
import Button from "@mui/material/Button";
import CircularProgress from "@mui/material/CircularProgress";
import Divider from "@mui/material/Divider";
import Stack from "@mui/material/Stack";
import Tooltip from "@mui/material/Tooltip";
import ToggleButton from "@mui/material/ToggleButton";
import ToggleButtonGroup from "@mui/material/ToggleButtonGroup";
import Typography from "@mui/material/Typography";
import DownloadIcon from "@mui/icons-material/Download";
import { useQuery } from "@tanstack/react-query";
import { DialogWindow } from "../DialogWindow";
import type { DiscoveryQueryItem } from "../../features/session/types";
import { DISCOVERY_SCORE_MODE_OPTIONS } from "../../features/discovery/types";
import { getDiscoveryQueryResult } from "../../features/discovery/api";
import { AlignmentGrid, ResultRow, DownloadSvgButton, sanitizeFilenamePart, type AlignmentGridRow } from "./AlignmentGrid";
import { WorkspaceCompare } from "./WorkspaceCompare";
import { downloadJson } from "./downloadJson";

type ViewMode = "pairwise" | "msa" | "compare";

type DialogViewDiscoveryQueryProps = {
  item: DiscoveryQueryItem;
  // The item's own session, or null when viewing an already-self-contained uploaded
  // result file (see the module doc comment above).
  sessionId: string | null;
  open: boolean;
  onClose: () => void;
};

function disabledReason(flagLabel: string): string {
  return `Not computed for this query. Submit a new query with "${flagLabel}" enabled.`;
}

export const DialogViewDiscoveryQuery: React.FC<DialogViewDiscoveryQueryProps> = ({ item, sessionId, open, onClose }) => {
  const [resultsView, setResultsView] = React.useState<ViewMode>("pairwise");
  const [selectedForMsa, setSelectedForMsa] = React.useState<Set<string>>(new Set());

  const resultQuery = useQuery({
    queryKey: ["getDiscoveryQueryResult", sessionId, item.id],
    queryFn: ({ signal }) => getDiscoveryQueryResult(sessionId as string, item.id, signal),
    enabled: open && sessionId !== null && item.status === "done",
    // Once computed, a query's payload never changes -- no need to ever refetch it.
    staleTime: Infinity,
  });

  const payload = sessionId === null ? item.payload ?? null : resultQuery.data ?? null;
  const results = React.useMemo(() => payload?.results ?? [], [payload]);

  // Reset the view mode each time the dialog opens for a (possibly different) item.
  React.useEffect(() => {
    if (!open) return;
    setResultsView("pairwise");
  }, [open, item.id]);

  // Default selection to "everyone in", same as the old single-shot flow's onSuccess
  // did -- keyed on `results` itself (not just `open`/`item.id`), since for a live
  // item `results` starts empty and only becomes real once resultQuery resolves;
  // payload never changes once fetched (staleTime: Infinity), so this only re-runs
  // when it actually has something new to select.
  React.useEffect(() => {
    if (!open) return;
    setSelectedForMsa(new Set(results.map((r) => r.entryId)));
  }, [open, results]);

  const toggleSelectedForMsa = React.useCallback((entryId: string) => {
    setSelectedForMsa((prev) => {
      const next = new Set(prev);
      if (next.has(entryId)) next.delete(entryId);
      else next.add(entryId);
      return next;
    });
  }, []);

  const selectedResults = React.useMemo(
    () => results.filter((r) => selectedForMsa.has(r.entryId)),
    [results, selectedForMsa]
  );

  // Precomputed MSA covers every top-X result, not just the current selection --
  // filtering the already-aligned rows down to the selection is free (same
  // coordinate system throughout), unlike the old flow which had to re-request a
  // fresh MSA from the server every time the selection changed.
  const msaRows: AlignmentGridRow[] = React.useMemo(() => {
    if (!payload?.msaRows) return [];
    const resultsById = new Map(results.map((r) => [r.entryId, r]));
    return payload.msaRows
      .filter((row) => row.id === "query" || selectedForMsa.has(row.id))
      .map((row) => {
        if (row.id === "query") {
          return { id: "query", label: "Query", sequence: row.alignedSequence };
        }
        const result = resultsById.get(row.id);
        return {
          id: row.id,
          label: result?.name ?? row.id,
          score: result ? `${result.normalizedAlignmentScorePct.toFixed(1)}%` : undefined,
          matchStrengths: row.similarityToQuery,
          sequence: row.alignedSequence,
        };
      });
  }, [payload, results, selectedForMsa]);

  const queryOriginSmiles = item.settings.queryOriginSmiles;
  const scoreModeLabel = DISCOVERY_SCORE_MODE_OPTIONS.find((o) => o.value === item.settings.scoreMode)?.label;

  return (
    <DialogWindow
      open={open}
      onClose={onClose}
      title={item.name}
      dividers
      maxWidth="lg"
      actions={[
        {
          key: "download",
          label: "Download result JSON",
          variant: "contained",
          color: "primary",
          disabled: item.status !== "done" || !payload,
          startIcon: <DownloadIcon fontSize="small" />,
          // Downloads the fully resolved item (payload included) -- not the raw `item`
          // prop, whose payload is stripped for a live session item (see module doc
          // comment). This is what makes the downloaded file self-contained enough to
          // re-upload and view later with no backend involved.
          onClick: () => downloadJson({ ...item, payload }, `discovery_query_${sanitizeFilenamePart(item.name)}.json`),
        },
        { key: "close", label: "Close", variant: "text", color: "inherit", onClick: onClose },
      ]}
    >
      {item.status === "queued" && (
        <Typography variant="body2" color="text.secondary">
          Waiting to run...
        </Typography>
      )}

      {item.status === "processing" && <CircularProgress size={24} />}

      {item.status === "error" && <Alert severity="error">{item.errorMessage || "Query failed."}</Alert>}

      {item.status === "done" && resultQuery.isLoading && <CircularProgress size={24} />}

      {item.status === "done" && resultQuery.error && (
        <Alert severity="error">{(resultQuery.error as Error).message || "Failed to load query results."}</Alert>
      )}

      {item.status === "done" && payload && (
        <Box>
          <Typography variant="body2" color="text.secondary" sx={{ mb: 1.5 }}>
            {payload.candidatesConsidered} candidate(s) considered
            {payload.candidatesSkipped > 0 ? `, ${payload.candidatesSkipped} skipped (unrecognized tokens)` : ""}
            . Showing top {payload.results.length}
            {scoreModeLabel ? ` (${scoreModeLabel.toLowerCase()})` : ""}.
          </Typography>

          {payload.results.length === 0 ? (
            <Typography variant="body2" color="text.secondary">
              No results.
            </Typography>
          ) : (
            <>
              <Stack direction="row" alignItems="center" sx={{ mb: 1, justifyContent: "space-between" }}>
                <Stack direction="row" spacing={1.5} alignItems="center" sx={{ mb: 2 }}>
                  <ToggleButtonGroup
                    size="small"
                    exclusive
                    value={resultsView}
                    onChange={(_, value) => value && setResultsView(value)}
                  >
                    <ToggleButton value="pairwise">Pairwise</ToggleButton>

                    <Tooltip title={item.flags.computeMsa ? "" : disabledReason("Compute MSA")} arrow>
                      <span>
                        <ToggleButton value="msa" disabled={!item.flags.computeMsa}>
                          Multiple sequence alignment
                        </ToggleButton>
                      </span>
                    </Tooltip>

                    <Tooltip title={item.flags.computeCompare ? "" : disabledReason("Compute compound comparison")} arrow>
                      <span>
                        <ToggleButton value="compare" disabled={!item.flags.computeCompare}>
                          Compare
                        </ToggleButton>
                      </span>
                    </Tooltip>
                  </ToggleButtonGroup>

                  <Typography variant="caption" color="text.secondary">
                    {selectedForMsa.size} of {payload.results.length} selected
                  </Typography>
                </Stack>

                <Stack direction="row" spacing={1.5} alignItems="center">
                  <Button size="small" onClick={() => setSelectedForMsa(new Set(payload.results.map((r) => r.entryId)))}>
                    Select all
                  </Button>
                  <Button size="small" onClick={() => setSelectedForMsa(new Set())}>
                    Clear
                  </Button>
                </Stack>
              </Stack>

              {resultsView !== "compare" && (
                <Typography variant="body2" color="text.secondary" sx={{ display: "block", mb: 1.5 }}>
                  Cell shading the alignment reflects each motif's structural similarity to the query's aligned motif in
                  that column. Darker means a closer match, so a sequence's weak spots stand out at a glance.
                </Typography>
              )}

              {resultsView == "compare" && (
                <Typography variant="body2" color="text.secondary" sx={{ display: "block", mb: 1.5 }}>
                  Each violin shows the distribution of one metric across the selected results, comparing every candidate back
                  to the query. Pick which metrics to include below.
                </Typography>
              )}

              <Divider sx={{ mb: 1 }} />

              {resultsView === "pairwise" ? (
                payload.results.map((result, idx) => (
                  <ResultRow
                    key={result.entryId}
                    result={result}
                    rank={idx + 1}
                    selectedForMsa={selectedForMsa.has(result.entryId)}
                    onToggleSelectedForMsa={() => toggleSelectedForMsa(result.entryId)}
                  />
                ))
              ) : resultsView === "msa" ? (
                selectedResults.length === 0 ? (
                  <Typography variant="body2" color="text.secondary">
                    Select at least one result above to view the multiple sequence alignment.
                  </Typography>
                ) : (
                  <Box sx={{ py: 1 }}>
                    <Box sx={{ mt: 1, display: "flex", justifyContent: "flex-end" }}>
                      <DownloadSvgButton rows={msaRows} filename="msa_alignment.svg" />
                    </Box>
                    <AlignmentGrid rows={msaRows} />
                  </Box>
                )
              ) : selectedResults.length === 0 ? (
                <Typography variant="body2" color="text.secondary">
                  Select at least one result above to compare compounds.
                </Typography>
              ) : (
                <Box sx={{ py: 1 }}>
                  <WorkspaceCompare
                    results={selectedResults}
                    queryOriginSmiles={queryOriginSmiles}
                    precomputedTanimoto={
                      payload.compareValues && payload.compareRadius !== undefined && payload.compareNBits !== undefined
                        ? { radius: payload.compareRadius, nBits: payload.compareNBits, values: payload.compareValues }
                        : null
                    }
                  />
                </Box>
              )}
            </>
          )}
        </Box>
      )}
    </DialogWindow>
  );
};
