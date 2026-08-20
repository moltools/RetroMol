import React from "react";
import Alert from "@mui/material/Alert";
import Box from "@mui/material/Box";
import Button from "@mui/material/Button";
import Card from "@mui/material/Card";
import CardContent from "@mui/material/CardContent";
import Chip from "@mui/material/Chip";
import CircularProgress from "@mui/material/CircularProgress";
import DownloadIcon from "@mui/icons-material/Download";
import MenuItem from "@mui/material/MenuItem";
import Stack from "@mui/material/Stack";
import Table from "@mui/material/Table";
import TableBody from "@mui/material/TableBody";
import TableCell from "@mui/material/TableCell";
import TableContainer from "@mui/material/TableContainer";
import TableHead from "@mui/material/TableHead";
import TablePagination from "@mui/material/TablePagination";
import TableRow from "@mui/material/TableRow";
import TextField from "@mui/material/TextField";
import Tooltip from "@mui/material/Tooltip";
import Typography from "@mui/material/Typography";
import { useNotifications } from "../../../NotificationProvider";
import { getAnnotationTerms, getBrowseEntries, browseEntriesExportUrl, MAX_BROWSE_ENTRIES } from "../../../../features/browse/api";
import { groupSourcesByDatabase } from "../../../../features/sources";
import type { AnnotationTerm, BrowseEntry } from "../../../../features/browse/types";
import type { EntryType } from "../../../../features/enrichment/types";

const ALL_CATEGORY = "all";
const ALL_TERM = "all";

export const WorkspaceBrowse: React.FC = () => {
  const { pushNotification } = useNotifications();

  const [entryType, setEntryType] = React.useState<EntryType | "all">("all");
  const [terms, setTerms] = React.useState<AnnotationTerm[]>([]);
  const [category, setCategory] = React.useState<string>(ALL_CATEGORY);
  const [termId, setTermId] = React.useState<string>(ALL_TERM);

  const [entries, setEntries] = React.useState<BrowseEntry[]>([]);
  const [totalCount, setTotalCount] = React.useState(0);
  const [truncated, setTruncated] = React.useState(false);
  const [hasSearched, setHasSearched] = React.useState(false);
  const [loading, setLoading] = React.useState(false);
  const [error, setError] = React.useState<string | null>(null);

  const [page, setPage] = React.useState(0);
  const [rowsPerPage, setRowsPerPage] = React.useState(25);

  React.useEffect(() => {
    const controller = new AbortController();
    getAnnotationTerms(undefined, controller.signal)
      .then((resp) => setTerms(resp.terms))
      .catch((err) => {
        if (controller.signal.aborted) return;
        pushNotification(`Failed to load annotation terms: ${err instanceof Error ? err.message : String(err)}`, "error");
      });
    return () => controller.abort();
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, []);

  const handleBrowse = async (event?: React.FormEvent) => {
    event?.preventDefault();

    setLoading(true);
    setError(null);
    try {
      const resp = await getBrowseEntries(entryType, termId === ALL_TERM ? null : termId);
      setEntries(resp.entries);
      setTotalCount(resp.totalCount);
      setTruncated(resp.truncated);
      setPage(0);
      setHasSearched(true);
    } catch (err) {
      setError(err instanceof Error ? err.message : String(err));
    } finally {
      setLoading(false);
    }
  };

  const categories = React.useMemo(
    () => Array.from(new Set(terms.map((t) => t.category))).sort(),
    [terms]
  );

  const termOptions = React.useMemo(
    () => terms.filter((t) => category === ALL_CATEGORY || t.category === category),
    [terms, category]
  );

  const handleCategoryChange = (value: string) => {
    setCategory(value);
    setTermId(ALL_TERM); // term choices depend on category, so reset
  };

  const paginatedEntries = entries.slice(page * rowsPerPage, page * rowsPerPage + rowsPerPage);
  const exportUrl = browseEntriesExportUrl(entryType, termId === ALL_TERM ? null : termId);

  return (
    <Box sx={{ width: "100%", mx: "auto", display: "flex", flexDirection: "column", gap: "16px" }}>
      <Card variant="outlined">
        <CardContent>
          <Typography component="h1" variant="subtitle1">
            Browse annotations
          </Typography>
          <Typography variant="body1" sx={{ mb: 2 }}>
            Browse database entries and their annotations, and download the filtered set as a TSV file.
          </Typography>

          <Box component="form" onSubmit={handleBrowse}>
            <Stack direction={{ xs: "column", sm: "row" }} spacing={1}>
              <TextField
                select
                label="Type"
                size="small"
                sx={{ minWidth: 160 }}
                value={entryType}
                onChange={(e) => setEntryType(e.target.value as EntryType | "all")}
              >
                <MenuItem value="all">All</MenuItem>
                <MenuItem value="compound">Compounds</MenuItem>
                <MenuItem value="bgc">Gene clusters (BGCs)</MenuItem>
              </TextField>

              <TextField
                select
                label="Annotation category"
                size="small"
                sx={{ minWidth: 200 }}
                value={category}
                onChange={(e) => handleCategoryChange(e.target.value)}
              >
                <MenuItem value={ALL_CATEGORY}>All categories</MenuItem>
                {categories.map((c) => (
                  <MenuItem key={c} value={c}>{c}</MenuItem>
                ))}
              </TextField>

              <TextField
                select
                label="Annotation label"
                size="small"
                sx={{ minWidth: 220 }}
                value={termId}
                onChange={(e) => setTermId(e.target.value)}
              >
                <MenuItem value={ALL_TERM}>Any</MenuItem>
                {termOptions.map((t) => (
                  <MenuItem key={t.id} value={t.id}>
                    {t.rank ? `${t.rank}: ${t.label}` : t.label}
                  </MenuItem>
                ))}
              </TextField>
            </Stack>

            <Stack direction="row" spacing={2} alignItems="center" sx={{ mt: 2 }}>
              <Button type="submit" variant="contained" disabled={loading}>
                {loading ? <CircularProgress size={20} /> : "Browse"}
              </Button>

              <Tooltip
                title={
                  !hasSearched
                    ? "Browse first to enable download"
                    : truncated
                    ? `Downloads only the first ${MAX_BROWSE_ENTRIES.toLocaleString()} matching rows (of ${totalCount.toLocaleString()}) -- narrow your filters to get the rest`
                    : "Download all matching rows as a TSV file"
                }
                arrow
              >
                <span>
                  <Button
                    variant="outlined"
                    startIcon={<DownloadIcon />}
                    component="a"
                    href={exportUrl}
                    disabled={!hasSearched}
                  >
                    Download TSV
                  </Button>
                </span>
              </Tooltip>
            </Stack>
          </Box>
        </CardContent>
      </Card>

      {error && (
        <Alert severity="warning" variant="outlined">
          Couldn't load entries: {error}
        </Alert>
      )}

      {hasSearched && truncated && !loading && !error && (
        <Alert severity="info" variant="outlined">
          {totalCount.toLocaleString()} entries match these filters. Showing (and downloading) only the
          first {MAX_BROWSE_ENTRIES.toLocaleString()}. Narrow the filters to see/export the rest.
        </Alert>
      )}

      {hasSearched && !loading && !error && (
        <Card variant="outlined">
          <CardContent>
            <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
              {truncated
                ? `${entries.length.toLocaleString()} of ${totalCount.toLocaleString()} matching entries shown`
                : `${entries.length.toLocaleString()} matching entries`}
            </Typography>

            <TableContainer sx={{ maxHeight: 520 }}>
              <Table size="small" stickyHeader>
                <TableHead>
                  <TableRow>
                    <TableCell>Name</TableCell>
                    <TableCell>Type</TableCell>
                    <TableCell>SMILES / InChIKey</TableCell>
                    <TableCell>Sources</TableCell>
                    <TableCell>Phylogeny</TableCell>
                    <TableCell>Chemical class</TableCell>
                  </TableRow>
                </TableHead>
                <TableBody>
                  {paginatedEntries.map((entry) => (
                    <TableRow key={entry.id} hover>
                      <TableCell>{entry.name}</TableCell>
                      <TableCell>{entry.type === "bgc" ? "BGC" : "Compound"}</TableCell>
                      <TableCell sx={{ maxWidth: 220, overflow: "hidden", textOverflow: "ellipsis", whiteSpace: "nowrap" }}>
                        {entry.type === "compound" ? entry.smiles ?? "-" : entry.id}
                      </TableCell>
                      <TableCell>
                        {groupSourcesByDatabase(entry.sources).map((g) => (
                          <Tooltip
                            key={g.databaseName}
                            title={g.items.map((s) => s.name).join(", ")}
                          >
                            <Chip
                              label={g.count > 1 ? `${g.databaseName} ×${g.count}` : g.databaseName}
                              size="small"
                              sx={{ mr: 0.5 }}
                            />
                          </Tooltip>
                        ))}
                      </TableCell>
                      <TableCell>
                        {[entry.phylogenyType, entry.genus, entry.species].filter(Boolean).join(" › ") || "-"}
                      </TableCell>
                      <TableCell>
                        {entry.chemicalClasses.length > 0
                          ? entry.chemicalClasses.map((c) => (
                              <Chip key={c} label={c} size="small" sx={{ mr: 0.5 }} />
                            ))
                          : "-"}
                      </TableCell>
                    </TableRow>
                  ))}
                </TableBody>
              </Table>
            </TableContainer>

            <TablePagination
              component="div"
              count={entries.length}
              page={page}
              onPageChange={(_, newPage) => setPage(newPage)}
              rowsPerPage={rowsPerPage}
              onRowsPerPageChange={(e) => {
                setRowsPerPage(parseInt(e.target.value, 10));
                setPage(0);
              }}
              rowsPerPageOptions={[10, 25, 50, 100]}
            />
          </CardContent>
        </Card>
      )}
    </Box>
  );
};
