import React from "react";
import Alert from "@mui/material/Alert";
import Box from "@mui/material/Box";
import Card from "@mui/material/Card";
import CardContent from "@mui/material/CardContent";
import Chip from "@mui/material/Chip";
import CircularProgress from "@mui/material/CircularProgress";
import Collapse from "@mui/material/Collapse";
import Divider from "@mui/material/Divider";
import Stack from "@mui/material/Stack";
import TextField from "@mui/material/TextField";
import ToggleButton from "@mui/material/ToggleButton";
import ToggleButtonGroup from "@mui/material/ToggleButtonGroup";
import Typography from "@mui/material/Typography";
import ExpandMoreIcon from "@mui/icons-material/ExpandMore";
import ExpandLessIcon from "@mui/icons-material/ExpandLess";
import { useQuery } from "@tanstack/react-query";
import { fetchRuleSet } from "../../features/rules/api";
import type { MatchingRule, ReactionRule } from "../../features/rules/types";
import { CopyIconButton } from "../CopyIconButton";
import { MinimalIconButton } from "../MinimalIconButton";
import { MotifName } from "../MotifName";
import { horizontalScrollSx } from "../../theme/scrollbarSx";
import SmilesDrawerContainer from "../SmilesDrawerContainer.js";
import { DrawingAttribution } from "../DrawingAttribution";
import { ReactionScheme } from "./ReactionScheme";

const STRUCTURE_SIZE = 130;

type RuleKind = "matching" | "reaction";

function matchesQuery(haystack: string[], query: string): boolean {
  const q = query.trim().toLowerCase();
  if (!q) return true;
  return haystack.some((h) => h.toLowerCase().includes(q));
}

function PropsList({ props }: { props: Record<string, unknown> }) {
  const entries = Object.entries(props);
  if (entries.length === 0) return null;

  return (
    <Stack spacing={0.5} sx={{ mt: 1 }}>
      <Typography variant="caption" color="text.secondary">
        Properties
      </Typography>
      <Stack direction="row" spacing={1} useFlexGap flexWrap="wrap">
        {entries.map(([key, value]) => (
          <Chip
            key={key}
            size="small"
            variant="outlined"
            label={`${key}: ${Array.isArray(value) ? value.join(", ") : String(value)}`}
          />
        ))}
      </Stack>
    </Stack>
  );
}

function MatchingRuleRow({ rule }: { rule: MatchingRule }) {
  const [expanded, setExpanded] = React.useState(false);

  return (
    <Box sx={{ borderBottom: "1px solid", borderColor: "divider", py: 1 }}>
      <Stack
        direction="row"
        alignItems="center"
        spacing={1.5}
        sx={{ cursor: "pointer" }}
        onClick={() => setExpanded((e) => !e)}
      >
        <Box sx={{ flex: 1, minWidth: 0 }}>
          <Typography variant="body2" fontWeight={500} noWrap>
            <MotifName name={rule.name} />
          </Typography>
        </Box>

        {rule.terminal && <Chip label="Terminal" size="small" sx={{ fontSize: "0.7rem" }} />}
        {rule.stereochemistry && <Chip label="Stereochemistry" size="small" color="info" sx={{ fontSize: "0.7rem", color: "white" }} />}
        {rule.pseudonyms.length > 0 && (
          <Typography variant="caption" color="text.secondary" sx={{ minWidth: 0, maxWidth: 260 }} noWrap>
            {rule.pseudonyms.map((p, i) => (
              <React.Fragment key={p}>
                {i > 0 && ", "}
                <MotifName name={p} />
              </React.Fragment>
            ))}
          </Typography>
        )}

        <MinimalIconButton size="small" tabIndex={-1}>
          {expanded ? <ExpandLessIcon fontSize="small" /> : <ExpandMoreIcon fontSize="small" />}
        </MinimalIconButton>
      </Stack>

      <Collapse in={expanded} unmountOnExit>
        <Stack direction="row" spacing={2} sx={{ mt: 1, pl: 1 }}>
          <Box>
            <SmilesDrawerContainer
              identifier={`matching-rule-${rule.id}`}
              smiles={rule.displaySmiles || rule.smiles}
              size={STRUCTURE_SIZE}
            />
            <DrawingAttribution library="smiles-drawer" />
          </Box>

          <Stack spacing={1} sx={{ flex: 1, minWidth: 0 }}>
            <Box>
              <Stack direction="row" alignItems="center" spacing={0.5}>
                <Typography variant="caption" color="text.secondary">
                  SMILES
                </Typography>
                <CopyIconButton text={rule.smiles} label="SMILES" />
              </Stack>
              <Typography variant="body2" sx={{ fontFamily: "monospace", wordBreak: "break-all" }}>
                {rule.smiles}
              </Typography>
            </Box>

            {rule.displaySmiles && (
              <Box>
                <Stack direction="row" alignItems="center" spacing={0.5}>
                  <Typography variant="caption" color="text.secondary">
                    Display SMILES
                  </Typography>
                  <CopyIconButton text={rule.displaySmiles} label="display SMILES" />
                </Stack>
                <Typography variant="body2" sx={{ fontFamily: "monospace", wordBreak: "break-all" }}>
                  {rule.displaySmiles}
                </Typography>
              </Box>
            )}

            {rule.pseudonyms.length > 0 && (
              <Box>
                <Typography variant="caption" color="text.secondary" sx={{ display: "block" }}>
                  Pseudonyms
                </Typography>
                <Stack direction="row" spacing={1} useFlexGap flexWrap="wrap">
                  {rule.pseudonyms.map((p) => (
                    <Chip key={p} size="small" label={<MotifName name={p} />} />
                  ))}
                </Stack>
              </Box>
            )}

            <PropsList props={rule.props} />
          </Stack>
        </Stack>
      </Collapse>
    </Box>
  );
}

function ReactionRuleRow({ rule }: { rule: ReactionRule }) {
  const [expanded, setExpanded] = React.useState(false);

  return (
    <Box sx={{ borderBottom: "1px solid", borderColor: "divider", py: 1 }}>
      <Stack
        direction="row"
        alignItems="center"
        spacing={1.5}
        sx={{ cursor: "pointer" }}
        onClick={() => setExpanded((e) => !e)}
      >
        <Box sx={{ flex: 1, minWidth: 0 }}>
          <Typography variant="body2" fontWeight={500} noWrap>
            {rule.name}
          </Typography>
        </Box>

        {rule.allowedInBulk && <Chip label="Bulk-allowed" size="small" sx={{ fontSize: "0.7rem" }} />}

        <MinimalIconButton size="small" tabIndex={-1}>
          {expanded ? <ExpandLessIcon fontSize="small" /> : <ExpandMoreIcon fontSize="small" />}
        </MinimalIconButton>
      </Stack>

      <Collapse in={expanded} unmountOnExit>
        <Stack spacing={1.5} sx={{ mt: 1, pl: 1 }}>
          <Box>
            <Stack direction="row" alignItems="center" spacing={0.8}>
              <Typography variant="caption" color="text.secondary">
                Reaction SMARTS
              </Typography>
              <CopyIconButton text={rule.smarts} label="reaction SMARTS" />
            </Stack>
            <Typography variant="body2" sx={{ fontFamily: "monospace", wordBreak: "break-all" }}>
              {rule.smarts}
            </Typography>
          </Box>

          <Box>
            <Typography variant="caption" color="text.secondary" sx={{ display: "block", mb: 0.5 }}>
              Reaction scheme
            </Typography>
            <Box sx={{ ...horizontalScrollSx, pb: 1 }}>
              <ReactionScheme id={rule.id} smarts={rule.smarts} />
            </Box>
            <DrawingAttribution library="rdkit" />
          </Box>

          <PropsList props={rule.props} />
        </Stack>
      </Collapse>
    </Box>
  );
}

export const WorkspaceRules: React.FC = () => {
  const [ruleKind, setRuleKind] = React.useState<RuleKind>("matching");
  const [search, setSearch] = React.useState("");

  const ruleSetQuery = useQuery({
    queryKey: ["ruleSet"],
    queryFn: ({ signal }) => fetchRuleSet(signal),
    staleTime: Infinity,
    gcTime: Infinity,
  });

  const filteredMatchingRules = React.useMemo(() => {
    const rules = ruleSetQuery.data?.matchingRules ?? [];
    if (!search.trim()) return rules;
    return rules.filter((r) => matchesQuery([r.name, r.smiles, r.displaySmiles ?? "", ...r.pseudonyms], search));
  }, [ruleSetQuery.data, search]);

  const filteredReactionRules = React.useMemo(() => {
    const rules = ruleSetQuery.data?.reactionRules ?? [];
    if (!search.trim()) return rules;
    return rules.filter((r) => matchesQuery([r.name, r.smarts], search));
  }, [ruleSetQuery.data, search]);

  return (
    <Box sx={{ width: "100%", mx: "auto", display: "flex", flexDirection: "column", gap: "16px" }}>
      <Card variant="outlined">
        <CardContent>
          <Typography component="h1" variant="subtitle1">
            Rules
          </Typography>
          <Typography variant="body1">
            Browse the matching rules (the named building-block motifs the current default rule set recognizes) and
            reaction rules (the bond-breaking transformations used to decompose a molecule into those motifs) that make
            up the current default rule set.
          </Typography>
        </CardContent>
      </Card>

      <Card variant="outlined">
        <CardContent>
          <Stack direction="row" alignItems="center" justifyContent="space-between" spacing={2} sx={{ mb: 2 }} flexWrap="wrap">
            <ToggleButtonGroup
              size="small"
              exclusive
              value={ruleKind}
              onChange={(_, value) => value && setRuleKind(value)}
            >
              <ToggleButton value="matching">
                Matching rules{ruleSetQuery.data ? ` (${ruleSetQuery.data.matchingRules.length})` : ""}
              </ToggleButton>
              <ToggleButton value="reaction">
                Reaction rules{ruleSetQuery.data ? ` (${ruleSetQuery.data.reactionRules.length})` : ""}
              </ToggleButton>
            </ToggleButtonGroup>

            <TextField
              size="small"
              placeholder={ruleKind === "matching" ? "Search by name, pseudonym, or SMILES" : "Search by name or SMARTS"}
              value={search}
              onChange={(e) => setSearch(e.target.value)}
              sx={{ minWidth: 295 }}
            />
          </Stack>

          <Divider sx={{ mb: 1 }} />

          {ruleSetQuery.isLoading && (
            <Box sx={{ display: "flex", justifyContent: "center", py: 4 }}>
              <CircularProgress size={24} />
            </Box>
          )}

          {ruleSetQuery.error && (
            <Alert severity="error">{(ruleSetQuery.error as Error).message || "Failed to load the rule set."}</Alert>
          )}

          {ruleSetQuery.data && (
            <>
              {ruleKind === "matching" ? (
                filteredMatchingRules.length === 0 ? (
                  <Typography variant="body2" color="text.secondary">
                    No matching rules found.
                  </Typography>
                ) : (
                  <Box sx={{ maxHeight: 640, overflowY: "auto" }}>
                    {filteredMatchingRules.map((rule) => (
                      <MatchingRuleRow key={rule.id} rule={rule} />
                    ))}
                  </Box>
                )
              ) : filteredReactionRules.length === 0 ? (
                <Typography variant="body2" color="text.secondary">
                  No reaction rules found.
                </Typography>
              ) : (
                <Box sx={{ maxHeight: 640, overflowY: "auto" }}>
                  {filteredReactionRules.map((rule) => (
                    <ReactionRuleRow key={rule.id} rule={rule} />
                  ))}
                </Box>
              )}
            </>
          )}
        </CardContent>
      </Card>
    </Box>
  );
};
