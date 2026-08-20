import React from "react";
import Alert from "@mui/material/Alert";
import Box from "@mui/material/Box";
import Button from "@mui/material/Button";
import Card from "@mui/material/Card";
import CardContent from "@mui/material/CardContent";
import Checkbox from "@mui/material/Checkbox";
import Chip from "@mui/material/Chip";
import CircularProgress from "@mui/material/CircularProgress";
import MenuItem from "@mui/material/MenuItem";
import Stack from "@mui/material/Stack";
import Table from "@mui/material/Table";
import TableBody from "@mui/material/TableBody";
import TableCell from "@mui/material/TableCell";
import TableContainer from "@mui/material/TableContainer";
import TableHead from "@mui/material/TableHead";
import TableRow from "@mui/material/TableRow";
import TextField from "@mui/material/TextField";
import Tooltip from "@mui/material/Tooltip";
import Typography from "@mui/material/Typography";
import { useNotifications } from "../../../NotificationProvider";
import { searchEntries, runEnrichmentAnalysis, MAX_ENTRY_SEARCH_RESULTS, MAX_ENRICHMENT_SELECTION } from "../../../../features/enrichment/api";
import { groupSourcesByDatabase } from "../../../../features/sources";
import type { EntryType, EnrichmentResult, SearchEntry } from "../../../../features/enrichment/types";

const Q_VALUE_SIGNIFICANT = 0.05;

export const WorkspaceEnrichment: React.FC = () => {
  const { pushNotification } = useNotifications();

  const [query, setQuery] = React.useState("");
  const [entryType, setEntryType] = React.useState<EntryType | "all">("all");
  const [searchResults, setSearchResults] = React.useState<SearchEntry[]>([]);
  const [searching, setSearching] = React.useState(false);
  const [searchError, setSearchError] = React.useState<string | null>(null);

  const [selectedIds, setSelectedIds] = React.useState<Set<string>>(new Set());

  const [running, setRunning] = React.useState(false);
  const [results, setResults] = React.useState<EnrichmentResult[] | null>(null);

  const handleSearch = async (event?: React.FormEvent) => {
    event?.preventDefault();
    if (!query.trim()) return;

    setSearching(true);
    setSearchError(null);
    try {
      const resp = await searchEntries(query.trim(), entryType);
      setSearchResults(resp.results);
    } catch (err) {
      setSearchError(err instanceof Error ? err.message : String(err));
    } finally {
      setSearching(false);
    }
  };

  const toggleSelected = (id: string) => {
    setSelectedIds((prev) => {
      const next = new Set(prev);
      if (next.has(id)) {
        next.delete(id);
      } else {
        if (next.size >= MAX_ENRICHMENT_SELECTION) {
          pushNotification(`You can select at most ${MAX_ENRICHMENT_SELECTION} entries.`, "warning");
          return prev;
        }
        next.add(id);
      }
      return next;
    });
  };

  const handleSelectAllResults = () => {
    setSelectedIds((prev) => {
      const next = new Set(prev);
      for (const entry of searchResults) {
        if (next.size >= MAX_ENRICHMENT_SELECTION) break;
        next.add(entry.id);
      }
      return next;
    });
  };

  const handleClearSelection = () => setSelectedIds(new Set());

  const handleRunEnrichment = async () => {
    if (selectedIds.size === 0) return;

    setRunning(true);
    setResults(null);
    try {
      const resp = await runEnrichmentAnalysis(Array.from(selectedIds));
      setResults(resp.results);
    } catch (err) {
      const msg = err instanceof Error ? err.message : String(err);
      pushNotification(`Enrichment analysis failed: ${msg}`, "error");
    } finally {
      setRunning(false);
    }
  };

  return (
    <Box sx={{ width: "100%", mx: "auto", display: "flex", flexDirection: "column", gap: "16px" }}>
      <Card variant="outlined">
        <CardContent>
          <Typography component="h1" variant="subtitle1">
            Enrichment analysis
          </Typography>
          <Typography variant="body1" sx={{ mb: 2 }}>
            Search the database, select up to {MAX_ENRICHMENT_SELECTION} entries, then test whether that
            selection is enriched (or depleted) for any annotation label compared to the rest of the database.
          </Typography>

          <Box component="form" onSubmit={handleSearch}>
            <Stack direction={{ xs: "column", sm: "row" }} spacing={1}>
              <TextField
                label="Search by name"
                size="small"
                fullWidth
                value={query}
                onChange={(e) => setQuery(e.target.value)}
              />
              <TextField
                select
                label="Type"
                size="small"
                sx={{ minWidth: 140 }}
                value={entryType}
                onChange={(e) => setEntryType(e.target.value as EntryType | "all")}
              >
                <MenuItem value="all">All</MenuItem>
                <MenuItem value="compound">Compounds</MenuItem>
                <MenuItem value="bgc">Gene clusters (BGCs)</MenuItem>
              </TextField>
            </Stack>

            <Stack direction="row" spacing={2} alignItems="center" sx={{ mt: 2 }}>
              <Button type="submit" variant="contained" disabled={searching || !query.trim()}>
                {searching ? <CircularProgress size={20} /> : "Search"}
              </Button>
            </Stack>
          </Box>

          {searchError && (
            <Alert severity="warning" variant="outlined" sx={{ mt: 2 }}>
              Search failed: {searchError}
            </Alert>
          )}
        </CardContent>
      </Card>

      {searchResults.length > 0 && (
        <Card variant="outlined">
          <CardContent>
            <Stack direction="row" justifyContent="space-between" alignItems="center" sx={{ mb: 1 }}>
              <Typography component="h2" variant="subtitle1">
                Results ({searchResults.length}{searchResults.length >= MAX_ENTRY_SEARCH_RESULTS ? "+" : ""})
              </Typography>
              <Stack direction="row" spacing={1}>
                <Button size="small" onClick={handleSelectAllResults}>Select all</Button>
                <Button size="small" onClick={handleClearSelection} disabled={selectedIds.size === 0}>
                  Clear selection
                </Button>
              </Stack>
            </Stack>

            <Typography variant="caption" color="text.secondary">
              {selectedIds.size} / {MAX_ENRICHMENT_SELECTION} selected
            </Typography>

            <TableContainer sx={{ maxHeight: 400, mt: 1 }}>
              <Table size="small" stickyHeader>
                <TableHead>
                  <TableRow>
                    <TableCell padding="checkbox" />
                    <TableCell>Name</TableCell>
                    <TableCell>Type</TableCell>
                    <TableCell>Sources</TableCell>
                  </TableRow>
                </TableHead>
                <TableBody>
                  {searchResults.map((entry) => (
                    <TableRow key={entry.id} hover onClick={() => toggleSelected(entry.id)} sx={{ cursor: "pointer" }}>
                      <TableCell padding="checkbox">
                        <Checkbox checked={selectedIds.has(entry.id)} onClick={(e) => e.stopPropagation()} onChange={() => toggleSelected(entry.id)} />
                      </TableCell>
                      <TableCell>{entry.name}</TableCell>
                      <TableCell>{entry.type === "bgc" ? "BGC" : "Compound"}</TableCell>
                      <TableCell>
                        {groupSourcesByDatabase(entry.sources).map((g) => (
                          <Tooltip key={g.databaseName} title={g.items.map((s) => s.name).join(", ")}>
                            <Chip
                              label={g.count > 1 ? `${g.databaseName} ×${g.count}` : g.databaseName}
                              size="small"
                              sx={{ mr: 0.5 }}
                            />
                          </Tooltip>
                        ))}
                      </TableCell>
                    </TableRow>
                  ))}
                </TableBody>
              </Table>
            </TableContainer>

            <Box sx={{ mt: 2 }}>
              <Button variant="contained" disabled={selectedIds.size === 0 || running} onClick={handleRunEnrichment}>
                {running ? <CircularProgress size={20} /> : `Run enrichment on ${selectedIds.size} selected`}
              </Button>
            </Box>
          </CardContent>
        </Card>
      )}

      {results && (
        <Card variant="outlined">
          <CardContent>
            <Typography component="h2" variant="subtitle1" sx={{ mb: 1 }}>
              Enrichment results
            </Typography>

            {results.length === 0 ? (
              <Typography variant="body2" color="text.secondary">
                No annotated terms were found on the selected entries.
              </Typography>
            ) : (
              <TableContainer>
                <Table size="small">
                  <TableHead>
                    <TableRow>
                      <TableCell>Category</TableCell>
                      <TableCell>Label</TableCell>
                      <TableCell align="right">Selected</TableCell>
                      <TableCell align="right">Background</TableCell>
                      <TableCell align="right">Fold</TableCell>
                      <TableCell>Direction</TableCell>
                      <TableCell align="right">p-value</TableCell>
                      <TableCell align="right">q-value</TableCell>
                    </TableRow>
                  </TableHead>
                  <TableBody>
                    {results.map((r) => (
                      <TableRow key={r.termId} sx={r.qValue < Q_VALUE_SIGNIFICANT ? { fontWeight: 600 } : undefined}>
                        <TableCell>{r.category ?? "-"}</TableCell>
                        <TableCell>
                          {r.label}
                          {r.qValue < Q_VALUE_SIGNIFICANT && (
                            <Chip label="significant" color="primary" size="small" sx={{ ml: 1 }} />
                          )}
                        </TableCell>
                        <TableCell align="right">{r.selectedWithTerm} / {r.selectedTotal}</TableCell>
                        <TableCell align="right">{r.backgroundWithTerm} / {r.backgroundTotal}</TableCell>
                        <TableCell align="right">{r.foldEnrichment != null ? r.foldEnrichment.toFixed(2) : "-"}</TableCell>
                        <TableCell>{r.direction}</TableCell>
                        <TableCell align="right">{r.pValue.toExponential(2)}</TableCell>
                        <TableCell align="right">{r.qValue.toExponential(2)}</TableCell>
                      </TableRow>
                    ))}
                  </TableBody>
                </Table>
              </TableContainer>
            )}
          </CardContent>
        </Card>
      )}
    </Box>
  );
};
