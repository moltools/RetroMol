// Enrichment analysis over a Discovery query's own nearest neighbors -- reuses the
// exact same /api/enrichmentAnalysis endpoint (Fisher's exact + Benjamini-Hochberg,
// see gui/src/server/routes/enrichment.py) the old standalone Enrichment tab used,
// just fed from the results dialog's own checkbox selection instead of a separate
// manual entry-search step.
import React from "react";
import Alert from "@mui/material/Alert";
import Box from "@mui/material/Box";
import Button from "@mui/material/Button";
import Chip from "@mui/material/Chip";
import CircularProgress from "@mui/material/CircularProgress";
import Stack from "@mui/material/Stack";
import MuiTooltip from "@mui/material/Tooltip";
import Typography from "@mui/material/Typography";
import CheckCircleIcon from "@mui/icons-material/CheckCircle";
import CancelIcon from "@mui/icons-material/Cancel";
import DownloadIcon from "@mui/icons-material/Download";
import RefreshIcon from "@mui/icons-material/Refresh";
import { DataGrid, type GridColDef } from "@mui/x-data-grid";
import { useQuery } from "@tanstack/react-query";
import { runEnrichmentAnalysis } from "../../features/enrichment/api";
import { Q_VALUE_SIGNIFICANT, type EnrichmentResult } from "../../features/enrichment/types";

function formatScientific(value: number): string {
  return value === 0 ? "0" : value.toExponential(2);
}

// annotation_terms.category/rank/label come straight from the pipeline in raw form
// (e.g. "chemical_class", "chebi_biological_role") -- humanized here for display only,
// same "touch only presentation, never the underlying value" approach as
// features/database/format.ts's toSentenceCase.
function humanizeAnnotationText(value: string): string {
  const spaced = value.replace(/_/g, " ");
  const sentenceCased = spaced.length > 0 ? spaced.charAt(0).toUpperCase() + spaced.slice(1) : spaced;
  return sentenceCased.replace(/\bchebi\b/gi, "ChEBI");
}

function downloadTsv(rows: EnrichmentResult[], filename: string): void {
  const headers = [
    "category", "rank", "label",
    "selectedWithTerm", "selectedTotal", "backgroundWithTerm", "backgroundTotal",
    "foldEnrichment", "direction", "pValue", "qValue", "significant",
  ];
  const lines = [headers.join("\t")];
  for (const r of rows) {
    lines.push(
      [
        r.category ?? "",
        r.rank ?? "",
        r.label,
        r.selectedWithTerm,
        r.selectedTotal,
        r.backgroundWithTerm,
        r.backgroundTotal,
        r.foldEnrichment ?? "",
        r.direction,
        r.pValue,
        r.qValue,
        r.qValue < Q_VALUE_SIGNIFICANT ? "yes" : "no",
      ].join("\t")
    );
  }

  const blob = new Blob([lines.join("\n")], { type: "text/tab-separated-values;charset=utf-8" });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = filename;
  document.body.appendChild(a);
  a.click();
  document.body.removeChild(a);
  URL.revokeObjectURL(url);
}

const columns: GridColDef<EnrichmentResult>[] = [
  {
    field: "category",
    headerName: "Category",
    width: 140,
    valueFormatter: (value: string | null) => (value ? humanizeAnnotationText(value) : "—"),
  },
  {
    field: "rank",
    headerName: "Rank",
    width: 140,
    valueFormatter: (value: string | null) => (value ? humanizeAnnotationText(value) : "—"),
  },
  {
    field: "label",
    headerName: "Label",
    width: 220,
    valueFormatter: (value: string) => humanizeAnnotationText(value),
  },
  {
    field: "selectedWithTerm",
    headerName: "Selected",
    width: 110,
    type: "number",
    renderCell: (params) => `${params.row.selectedWithTerm} / ${params.row.selectedTotal}`,
  },
  {
    field: "backgroundWithTerm",
    headerName: "Background",
    width: 120,
    type: "number",
    renderCell: (params) => `${params.row.backgroundWithTerm} / ${params.row.backgroundTotal}`,
  },
  {
    field: "foldEnrichment",
    headerName: "Fold",
    width: 100,
    type: "number",
    valueFormatter: (value: number | null) => (value === null ? "—" : value.toFixed(2)),
  },
  { field: "direction", headerName: "Direction", width: 110 },
  {
    field: "pValue",
    headerName: "p-value",
    width: 110,
    type: "number",
    valueFormatter: (value: number) => formatScientific(value),
  },
  {
    field: "qValue",
    headerName: "q-value",
    width: 110,
    type: "number",
    valueFormatter: (value: number) => formatScientific(value),
  },
  {
    field: "significant",
    headerName: "Significant",
    width: 120,
    type: "boolean",
    valueGetter: (_value, row) => row.qValue < Q_VALUE_SIGNIFICANT,
    renderCell: (params) =>
      params.value ? (
        <CheckCircleIcon fontSize="small" color="success" titleAccess="Significant" />
      ) : (
        <CancelIcon fontSize="small" color="error" titleAccess="Not significant" />
      ),
  },
];

// DataGrid's own "baseTooltip" slot type has no `arrow` prop (it's a stripped-down
// interface, not @mui/material's full TooltipProps -- see gridBaseSlots.d.ts), so
// arrow-style tooltips (sort/filter/menu hover hints) need a full slot override
// rather than a slotProps tweak.
function GridArrowTooltip(props: React.ComponentProps<typeof MuiTooltip>) {
  return <MuiTooltip arrow {...props} />;
}

export const DiscoveryEnrichmentView: React.FC<{ entryIds: string[] }> = ({ entryIds }) => {
  // Recomputes only on mount (first open) and on explicit "Recalculate" clicks --
  // never automatically as the user toggles checkboxes above, which would otherwise
  // fire a request on every single click while they're still adjusting the selection.
  const [queriedIds, setQueriedIds] = React.useState<string[]>(entryIds);
  const didInit = React.useRef(false);
  React.useEffect(() => {
    if (didInit.current) return;
    didInit.current = true;
    setQueriedIds(entryIds);
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, []);

  const selectionChanged = React.useMemo(() => {
    const a = [...entryIds].sort();
    const b = [...queriedIds].sort();
    return a.length !== b.length || a.some((id, i) => id !== b[i]);
  }, [entryIds, queriedIds]);

  const enrichmentQuery = useQuery({
    queryKey: ["discoveryEnrichmentAnalysis", queriedIds],
    queryFn: ({ signal }) => runEnrichmentAnalysis(queriedIds, signal),
    enabled: queriedIds.length > 0,
  });

  const rows = enrichmentQuery.data?.results ?? [];

  return (
    <Box>
      <Stack direction="row" alignItems="center" justifyContent="space-between" flexWrap="wrap" gap={1} sx={{ mb: 1.5 }}>
        <Stack direction="row" alignItems="center" spacing={1}>
          <Typography variant="body2" color="text.secondary">
            Enrichment of the {queriedIds.length} selected nearest neighbor{queriedIds.length === 1 ? "" : "s"} vs. the
            rest of the database, one two-sided Fisher's exact test per annotation term, Benjamini-Hochberg corrected.
          </Typography>
          {selectionChanged && (
            <Chip
              label="Selection changed"
              size="small"
              color="warning"
              variant="outlined"
              sx={{ flexShrink: 0 }}
            />
          )}
        </Stack>
        <Stack direction="row" spacing={1} sx={{ flexShrink: 0 }}>
          <Button
            size="small"
            variant={selectionChanged ? "contained" : "outlined"}
            startIcon={enrichmentQuery.isFetching ? <CircularProgress size={14} color="inherit" /> : <RefreshIcon fontSize="small" />}
            onClick={() => setQueriedIds(entryIds)}
            disabled={entryIds.length === 0 || enrichmentQuery.isFetching}
          >
            Recalculate
          </Button>
          <Button
            size="small"
            startIcon={<DownloadIcon fontSize="small" />}
            onClick={() => downloadTsv(rows, "enrichment_analysis.tsv")}
            disabled={rows.length === 0}
          >
            Download TSV
          </Button>
        </Stack>
      </Stack>

      {queriedIds.length === 0 && (
        <Typography variant="body2" color="text.secondary">
          Select at least one result above, then click Recalculate to run enrichment analysis.
        </Typography>
      )}

      {enrichmentQuery.error && (
        <Alert severity="error" sx={{ mb: 1.5 }}>
          {(enrichmentQuery.error as Error).message || "Failed to run enrichment analysis."}
        </Alert>
      )}

      {queriedIds.length > 0 && (
        <Box sx={{ height: 480, width: "100%" }}>
          <DataGrid
            rows={rows}
            columns={columns}
            getRowId={(row) => row.termId}
            loading={enrichmentQuery.isLoading}
            density="compact"
            disableRowSelectionOnClick
            showCellVerticalBorder
            showColumnVerticalBorder
            initialState={{
              sorting: { sortModel: [{ field: "qValue", sort: "asc" }] },
            }}
            pageSizeOptions={[25, 50, 100]}
            slots={{ baseTooltip: GridArrowTooltip }}
            sx={{
              // Explicit border rules rather than relying on the --DataGrid-rowBorderColor
              // CSS variable / showCellVerticalBorder alone -- confirmed live that the
              // variable covers vertical cell borders but not the horizontal row
              // separator in this DataGrid version, so it's set directly here instead.
              // "divider" resolves through the app's own theme (light/dark both), not a
              // hardcoded color.
              "& .MuiDataGrid-cell": {
                fontSize: "0.8125rem",
                borderRight: "1px solid",
                borderRightColor: "divider",
                borderBottom: "1px solid",
                borderBottomColor: "divider",
              },
              "& .MuiDataGrid-columnHeader": {
                borderRight: "1px solid",
                borderRightColor: "divider",
              },
              // Sort arrow / column-menu (three dots) icons -- bare icon, no button
              // chrome (padding, hover box, border) around them.
              "& .MuiDataGrid-iconButtonContainer .MuiIconButton-root, & .MuiDataGrid-menuIconButton": {
                padding: 0,
                border: "none",
                backgroundColor: "transparent",
                boxShadow: "none",
                "&:hover": { backgroundColor: "transparent" },
              },
            }}
          />
        </Box>
      )}
    </Box>
  );
};
