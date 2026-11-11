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
import MenuItem from "@mui/material/MenuItem";
import Select from "@mui/material/Select";
import Stack from "@mui/material/Stack";
import Tooltip from "@mui/material/Tooltip";
import Typography from "@mui/material/Typography";
import DownloadIcon from "@mui/icons-material/Download";
import { useQuery } from "@tanstack/react-query";
import { DialogWindow } from "../DialogWindow";
import type { Session, DiscoveryQueryItem } from "../../features/session/types";
import { DISCOVERY_SCORE_MODE_OPTIONS, type DiscoveryResult } from "../../features/discovery/types";
import { getDiscoveryQueryResult } from "../../features/discovery/api";
import { AlignmentGrid, ResultRow, DownloadSvgButton, sanitizeFilenamePart, type AlignmentGridRow } from "./AlignmentGrid";
import { WorkspaceCompare } from "./WorkspaceCompare";
import { DiscoveryEnrichmentView } from "./DiscoveryEnrichmentView";
import { downloadJson } from "./downloadJson";
import { useNotifications } from "../NotificationProvider";
import { MAX_ITEMS, importCompound, importClustersBatch } from "../../features/jobs/api";

type ViewMode = "pairwise" | "msa" | "compare" | "enrichment";

const VIEW_MODE_OPTIONS: {
  value: ViewMode;
  label: string;
  flag: keyof DiscoveryQueryItem["flags"] | null;
  // The "Compute for this query" checkbox's own label (see WorkspaceDiscovery.tsx),
  // for the disabled-option tooltip -- deliberately not just reusing `label` above,
  // since the view's display name ("Multiple sequence alignment") and the checkbox
  // that enables it ("Compute MSA") aren't worded the same.
  checkboxLabel: string | null;
}[] = [
  { value: "pairwise", label: "Pairwise", flag: null, checkboxLabel: null },
  { value: "msa", label: "Multiple sequence alignment", flag: "computeMsa", checkboxLabel: "Compute MSA" },
  { value: "compare", label: "Compare", flag: "computeCompare", checkboxLabel: "Compute compound comparison" },
  { value: "enrichment", label: "Enrichment", flag: "computeEnrichment", checkboxLabel: "Compute enrichment" },
];

// Same default the "Import BGCs" dialog offers -- a sent-back BGC result has no UI of
// its own to pick a threshold, so it's reparsed with the standard starting point.
const DEFAULT_PARAS_THRESHOLD = 0.1;

type DialogViewDiscoveryQueryProps = {
  item: DiscoveryQueryItem;
  // The item's own session, or null when viewing an already-self-contained uploaded
  // result file (see the module doc comment above).
  sessionId: string | null;
  // The live workspace session/setter -- always the *current* session regardless of
  // where the viewed query came from, since "send back to Uploads" adds to whatever
  // session is open right now, not necessarily the one the query originally ran in.
  session: Session;
  setSession: React.Dispatch<React.SetStateAction<Session | null>>;
  open: boolean;
  onClose: () => void;
};

function disabledReason(flagLabel: string): string {
  return `Not computed for this query. Submit a new query with "${flagLabel}" enabled.`;
}

export const DialogViewDiscoveryQuery: React.FC<DialogViewDiscoveryQueryProps> = ({
  item,
  sessionId,
  session,
  setSession,
  open,
  onClose,
}) => {
  const [resultsView, setResultsView] = React.useState<ViewMode>("pairwise");
  const [selectedForMsa, setSelectedForMsa] = React.useState<Set<string>>(new Set());
  const [sendingEntryIds, setSendingEntryIds] = React.useState<Set<string>>(new Set());
  const { pushNotification } = useNotifications();

  const setSessionSafe = React.useCallback(
    (updater: (prev: Session) => Session) => {
      setSession((prev) => (prev ? updater(prev) : prev));
    },
    [setSession]
  );

  const importDeps = React.useMemo(
    () => ({ setSession: setSessionSafe, pushNotification, sessionId: session.sessionId }),
    [setSessionSafe, pushNotification, session.sessionId]
  );

  const uploadItemCount = React.useMemo(
    () => session.items.filter((it) => it.kind === "compound" || it.kind === "cluster").length,
    [session.items]
  );
  const uploadSlotsFull = uploadItemCount >= MAX_ITEMS;

  const handleSendToUploads = React.useCallback(
    async (result: DiscoveryResult) => {
      if (!result.raw) return;

      setSendingEntryIds((prev) => new Set(prev).add(result.entryId));
      try {
        const sent =
          result.type === "compound"
            ? await importCompound(importDeps, { name: result.name, smiles: result.raw, matchStereochemistry: false })
            : (await importClustersBatch(importDeps, [
                { name: result.name, fileContent: result.raw, parasThreshold: DEFAULT_PARAS_THRESHOLD },
              ]))[0] ?? null;

        if (sent) {
          pushNotification(`Sent "${result.name}" to Uploads for reparsing.`, "success");
        }
      } finally {
        setSendingEntryIds((prev) => {
          const next = new Set(prev);
          next.delete(result.entryId);
          return next;
        });
      }
    },
    [importDeps, pushNotification]
  );

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
                  <Select
                    size="small"
                    value={resultsView}
                    onChange={(e) => setResultsView(e.target.value as ViewMode)}
                    sx={{ minWidth: 220 }}
                  >
                    {VIEW_MODE_OPTIONS.map((option) => {
                      const enabled = option.flag === null || item.flags[option.flag];
                      const menuItem = (
                        <MenuItem key={option.value} value={option.value} disabled={!enabled}>
                          {option.label}
                        </MenuItem>
                      );
                      if (enabled) return menuItem;
                      // A disabled MenuItem still needs a wrapping span for the tooltip to
                      // fire (disabled elements don't receive pointer events on their own).
                      return (
                        <Tooltip
                          key={option.value}
                          title={disabledReason(option.checkboxLabel ?? option.label)}
                          arrow
                          placement="right"
                        >
                          <span>{menuItem}</span>
                        </Tooltip>
                      );
                    })}
                  </Select>

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

              {(resultsView === "pairwise" || resultsView === "msa") && (
                <Typography variant="body2" color="text.secondary" sx={{ display: "block", mb: 1.5 }}>
                  Cell shading the alignment reflects each motif's structural similarity to the query's aligned motif in
                  that column. Darker means a closer match, so a sequence's weak spots stand out at a glance.
                </Typography>
              )}

              {resultsView === "compare" && (
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
                    onSendToUploads={handleSendToUploads}
                    uploadSlotsFull={uploadSlotsFull}
                    sendingToUploads={sendingEntryIds.has(result.entryId)}
                  />
                ))
              ) : resultsView === "msa" ? (
                selectedResults.length === 0 ? (
                  <Typography variant="body2" color="text.secondary">
                    Select at least one result above to view the multiple sequence alignment.
                  </Typography>
                ) : (
                  <Box sx={{ py: 1 }}>
                    <Box sx={{ mt: 1, mb: 1.5, display: "flex", justifyContent: "flex-end" }}>
                      <DownloadSvgButton rows={msaRows} filename="msa_alignment.svg" />
                    </Box>
                    <AlignmentGrid rows={msaRows} />
                  </Box>
                )
              ) : resultsView === "enrichment" ? (
                <Box sx={{ py: 1 }}>
                  <DiscoveryEnrichmentView entryIds={Array.from(selectedForMsa)} />
                </Box>
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
