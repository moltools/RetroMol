import React from "react";
import Box from "@mui/material/Box";
import Card from "@mui/material/Card";
import CardContent from "@mui/material/CardContent";
import Typography from "@mui/material/Typography";
import MuiLink from "@mui/material/Link";
import Stack from "@mui/material/Stack";
import InputLabel from "@mui/material/InputLabel";
import FormControl from "@mui/material/FormControl";
import MenuItem from "@mui/material/MenuItem";
import Button from "@mui/material/Button";
import CircularProgress from '@mui/material/CircularProgress';
import Alert from "@mui/material/Alert";
import Checkbox from "@mui/material/Checkbox";
import FormGroup from "@mui/material/FormGroup";
import FormControlLabel from "@mui/material/FormControlLabel";
import { useTheme } from "@mui/material/styles";
import { useNotifications } from "../NotificationProvider";
import { Link as RouterLink } from "react-router-dom";
import { Session } from "../../features/session/types";
import { Select } from "@mui/material";
import { QueryResultView } from "./QueryResultView";
import { QueryResult } from "./query_types";

type WorkspaceDiscoveryProps = {
  session: Session;
  setSession: React.Dispatch<React.SetStateAction<Session | null>>;
};

export const WorkspaceDiscovery: React.FC<WorkspaceDiscoveryProps> = ({ session, setSession }) => {
  const theme = useTheme();
  const { pushNotification } = useNotifications();

  // Check if session has items
  const hasItems = session.items.length > 0;

  // Query state
  const [selectedItemId, setSelectedItemId] = React.useState<string>("");
  const [queryLoading, setQueryLoading] = React.useState(false);
  const [queryError, setQueryError] = React.useState<string | null>(null);
  const [queryResult, setQueryResult] = React.useState<QueryResult | null>(null);

  // Query settings
  const [queryAgainstCompounds, setQueryAgainstCompounds] = React.useState(true);
  const [queryAgainstClusters, setQueryAgainstClusters] = React.useState(false);
  const [queryAgainstUserUploads, setQueryAgainstUserUploads] = React.useState(false);

  // Wrap parent setter (Session | null) into the deps shape (Session-only functional updater)
  const setSessionSafe = React.useCallback(
    (updater: (prev: Session) => Session) => {
      setSession((prev) => (prev ? updater(prev) : prev));
    },
    [setSession]
  );

  // Memoized alert based on query state
  const alert = React.useMemo(() => {
    if (queryError) {
      return { severity: "error" as const, text: queryError };
    }
    if (queryResult) {
      return { severity: "success" as const, text: "Query complete! See results below." };
    }
    return {
      severity: "info" as const,
      text: "Select an item and click “Run query” to see results here.",
    };
  }, [queryError, queryResult]);

  // Post query
  async function queryItem(itemId: string): Promise<QueryResult> {
    const params = new URLSearchParams({
      sessionId: session.sessionId,
      itemId,
      queryAgainstUserUploads: String(queryAgainstUserUploads),
      queryAgainstCompounds: String(queryAgainstCompounds),
      queryAgainstClusters: String(queryAgainstClusters)
    });
    const res = await fetch(`/api/queryItem?${params.toString()}`);
    if (!res.ok) { throw new Error(`Query failed: ${res.status}`); };
    return await res.json();
  };

  // Handler to run query (dummy implementation)
  const handleRunQuery = async () => {
    if (!selectedItemId) return;

    setQueryLoading(true);
    setQueryError(null);
    setQueryResult(null);

    try {
      const result = await queryItem(selectedItemId);
      setQueryResult(result);
      pushNotification("Query completed successfully!", "success");
    } catch (error) {
      setQueryError("Failed to run query. Please try again.");
      pushNotification("Failed to run query.", "error");
    } finally {
      setQueryLoading(false);
    }
  };

  return (
    <Box
      sx={{
          width: "100%",
          mx: "auto",
          display: "flex",
          flexDirection: "column",
          gap: "16px"
      }}
    >
      <Card
        variant="outlined"
        sx={{
          display: "flex",
          flexDirection: "column",
          gap: "8px",
          flexGrow: 1,
        }}
      >
        <CardContent>
          <Typography component="h1" variant="subtitle1">
            Getting started
          </Typography>
          <Typography variant="body1">
            Here you can use uploaded items from the&nbsp;
            <MuiLink
              component={RouterLink}
              to="/dashboard/upload"
              underline="hover"
              color={(theme.vars || theme).palette.primary.main}
              sx={{ fontWeight: "500" }}
            >
              Upload tab
            </MuiLink>
            &nbsp;for cross-modal retrieval against the BioNexus database or other items in your workspace.
          </Typography>

          <Stack direction="column" spacing={2} sx={{ mt: 3 }}>
            <FormControl fullWidth size="small" disabled={!hasItems || queryLoading}>
              <InputLabel
                id="item-select-label"
                sx={{
                  backgroundColor: "background.paper",
                  px: 0.5,
                  transform:
                    // We adjust the height to vertically center with the Select input
                    "translate(14px, 12px) scale(1)",
                    "&.MuiInputLabel-shrink": { transform: "translate(14px, -9px) scale(0.75)" },
                }}
              >
                {!hasItems
                  ? "No items available to select"
                  : selectedItemId
                    ? "Item to use for querying"
                    : "Select an item to use for querying"}
              </InputLabel>
              <Select
                labelId="item-select-label"
                label="Item to use for querying"
                value={selectedItemId}
                onChange={(e) => setSelectedItemId(e.target.value)}
                disabled={!hasItems || queryLoading}
                MenuProps={{
                  PaperProps: {
                    sx: {
                      // Prevent "drag highlight" in the dropdown options
                      "& .MuiMenuItem-root": {
                        userSelect: "none",
                        borderRadius: 0,
                      },
                    },
                  },
                }}
                sx={{
                  // Prevent "drag highlight" on the selected value area
                  "& .MuiSelect-select": { userSelect: "none" },
                  // MUI adds a focus background to the select display: remove it!
                  "& .MuiSelect-select:focus": { backgroundColor: "transparent" },
                  // This usually removes the browser-style outline
                  "& .MuiSelect-select:focus-visible": { outline: "none" },
                  // Make just as high as Buttons
                  "&.MuiInputBase-root": { height: 44 },
                }}
              >
                {session.items.map((item) => (
                  <MenuItem key={item.id} value={item.id}>
                    {item.name}
                  </MenuItem>
                ))}
              </Select>
            </FormControl>

            <Box sx={{ display: "flex", justifyContent: "center" }}>
              <FormGroup row>
                {/*<FormControlLabel*/}
                {/*  control={*/}
                {/*    <Checkbox*/}
                {/*      checked={queryAgainstClusters}*/}
                {/*      disabled={*/}
                {/*        queryLoading ||*/}
                {/*        (queryAgainstClusters && !queryAgainstCompounds && !queryAgainstUserUploads)*/}
                {/*      }*/}
                {/*      onChange={(e) => setQueryAgainstClusters(e.target.checked)}*/}
                {/*    />*/}
                {/*  }*/}
                {/*  label="Query against clusters"*/}
                {/*/>*/}
                <FormControlLabel
                  control={
                    <Checkbox
                      checked={queryAgainstCompounds}
                      disabled={true}
                      // disabled={
                      //   queryLoading ||
                      //   (queryAgainstCompounds && !queryAgainstClusters && !queryAgainstUserUploads)
                      // }
                      onChange={(e) => setQueryAgainstCompounds(e.target.checked)}
                    />
                  }
                  label="Query against compounds"
                />
                {/*<FormControlLabel*/}
                {/*  control={*/}
                {/*    <Checkbox*/}
                {/*      checked={queryAgainstUserUploads}*/}
                {/*      disabled={*/}
                {/*        queryLoading ||*/}
                {/*        (queryAgainstUserUploads && !queryAgainstClusters && !queryAgainstCompounds)*/}
                {/*      }*/}
                {/*      onChange={(e) => setQueryAgainstUserUploads(e.target.checked)}*/}
                {/*    />*/}
                {/*  }*/}
                {/*  label="Query against user uploads"*/}
                {/*/>*/}
              </FormGroup>
            </Box>

            <Stack direction="row" spacing={1} alignItems="center">
              <Button variant="contained" onClick={handleRunQuery} disabled={queryLoading || !selectedItemId}>
                {queryLoading ? "Running query..." : "Run query"}
              </Button>
              {queryLoading && <CircularProgress size={24} />}
            </Stack>

            <Alert severity={alert.severity}>{alert.text}</Alert>
          </Stack>

        </CardContent>
      </Card>

      {queryResult && (
        <Card
          variant="outlined"
          sx={{
            display: "flex",
            flexDirection: "column",
            gap: "8px",
            flexGrow: 1,
          }}
        >
          <CardContent>
            <QueryResultView result={queryResult} />
          </CardContent>
        </Card>
      )}
    </Box>
  );
};