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
import { getAnnotationTerms, getBrowseEntries, browseEntriesExportUrl } from "../../../../features/browse/api";
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
  const [loading, setLoading] = React.useState(true);
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

  React.useEffect(() => {
    const controller = new AbortController();
    setLoading(true);
    setError(null);

    getBrowseEntries(entryType, termId === ALL_TERM ? null : termId, controller.signal)
      .then((resp) => {
        setEntries(resp.entries);
        setPage(0);
      })
      .catch((err) => {
        if (controller.signal.aborted) return;
        setError(err instanceof Error ? err.message : String(err));
      })
      .finally(() => {
        if (controller.signal.aborted) return;
        setLoading(false);
      });

    return () => controller.abort();
  }, [entryType, termId]);

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
            Browse database entries and their annotations, and download the filtered set as a TSV file
            (compound SMILES/InChIKey, or BGC id, alongside phylogeny and chemical class annotations).
          </Typography>

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

            <Tooltip title="Download the entries below (all matching rows, not just this page) as a TSV file">
              <span>
                <Button
                  variant="contained"
                  startIcon={<DownloadIcon />}
                  component="a"
                  href={exportUrl}
                  disabled={loading || entries.length === 0}
                >
                  Download TSV
                </Button>
              </span>
            </Tooltip>
          </Stack>
        </CardContent>
      </Card>

      {error && (
        <Alert severity="warning" variant="outlined">
          Couldn't load entries: {error}
        </Alert>
      )}

      {loading && (
        <Box sx={{ display: "flex", justifyContent: "center", py: 4 }}>
          <CircularProgress />
        </Box>
      )}

      {!loading && !error && (
        <Card variant="outlined">
          <CardContent>
            <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
              {entries.length.toLocaleString()} matching entries
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
                        {entry.sources.map((s) => (
                          <Chip key={s.databaseName + s.name} label={s.databaseName} size="small" sx={{ mr: 0.5 }} />
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
