import React from "react";
import Alert from "@mui/material/Alert";
import Box from "@mui/material/Box";
import Button from "@mui/material/Button";
import Checkbox from "@mui/material/Checkbox";
import Chip from "@mui/material/Chip";
import CircularProgress from "@mui/material/CircularProgress";
import Divider from "@mui/material/Divider";
import FormControlLabel from "@mui/material/FormControlLabel";
import MenuItem from "@mui/material/MenuItem";
import Popover from "@mui/material/Popover";
import Stack from "@mui/material/Stack";
import TextField from "@mui/material/TextField";
import Typography from "@mui/material/Typography";
import ArrowDropDownIcon from "@mui/icons-material/ArrowDropDown";
import DownloadIcon from "@mui/icons-material/Download";
import TuneIcon from "@mui/icons-material/Tune";
import { useQuery } from "@tanstack/react-query";
import { alpha, useTheme, type Theme } from "@mui/material/styles";
import { runDiscoveryCompareTanimoto } from "../../features/discovery/api";
import {
  MORGAN_NBITS_DEFAULT,
  MORGAN_NBITS_OPTIONS,
  MORGAN_RADIUS_DEFAULT,
  type DiscoveryResult,
} from "../../features/discovery/types";
import { exportChartSvg } from "./chartSvgExport";

type CompareMetricId = "fingerprintCosine" | "normalizedAlignment" | "tanimoto";

// The single source of truth for what metrics the picker offers -- adding a metric
// (sync or server-backed) means adding one entry here plus a branch where `metrics`
// is built below, never touching the picker UI itself. Order here also fixes each
// metric's categorical color slot (see chartColorFor) and its position in the plot.
type MetricDef = {
  id: CompareMetricId;
  label: string;
  description: string;
};

const METRIC_DEFS: MetricDef[] = [
  {
    id: "fingerprintCosine",
    label: "Fingerprint cosine",
    description: "Cosine similarity between per-monomer token fingerprints -- the same metric Discovery ranks results by.",
  },
  {
    id: "normalizedAlignment",
    label: "Normalized alignment",
    description: "Sequence alignment score, normalized to a percentage of the best possible score.",
  },
  {
    id: "tanimoto",
    label: "Tanimoto (structure)",
    description: "Real molecular structure similarity via Morgan/ECFP fingerprints -- only compounds with a known SMILES.",
  },
];

const METRIC_LABELS: Record<CompareMetricId, string> = Object.fromEntries(
  METRIC_DEFS.map((def) => [def.id, def.label])
) as Record<CompareMetricId, string>;

type CompareDataPoint = {
  id: string;
  label: string;
  value: number | null; // fraction in [0, 1]; null = no data for this candidate
};

type CompareMetricResult = {
  id: CompareMetricId;
  label: string;
  color: string;
  points: CompareDataPoint[];
  isLoading?: boolean;
  error?: string | null;
};

// ---- Violin geometry -------------------------------------------------------

const KDE_SAMPLES = 41; // sampled at 0, 0.025, ..., 1
const SINGLE_POINT_BANDWIDTH = 0.045;
const MIN_BANDWIDTH = 0.02;

function gaussianKdeAt(values: number[], bandwidth: number, x: number): number {
  const n = values.length;
  let sum = 0;
  for (const v of values) {
    const z = (x - v) / bandwidth;
    sum += Math.exp(-0.5 * z * z);
  }
  return sum / (n * bandwidth * Math.sqrt(2 * Math.PI));
}

// Samples a Gaussian KDE of `values` (each expected in [0, 1]) at KDE_SAMPLES evenly
// spaced points across [0, 1] -- the fixed domain every metric here shares, since
// they're all similarity/score fractions. Bandwidth follows Silverman's rule of
// thumb, floored so a tight cluster (or a single value) still renders a visible
// shape instead of collapsing to a needle.
function sampleKde(values: number[]): { xs: number[]; density: number[] } {
  const xs = Array.from({ length: KDE_SAMPLES }, (_, i) => i / (KDE_SAMPLES - 1));

  if (values.length === 0) return { xs, density: xs.map(() => 0) };

  let bandwidth: number;
  if (values.length === 1) {
    bandwidth = SINGLE_POINT_BANDWIDTH;
  } else {
    const mean = values.reduce((a, b) => a + b, 0) / values.length;
    const variance = values.reduce((a, b) => a + (b - mean) ** 2, 0) / (values.length - 1);
    const std = Math.sqrt(variance);
    bandwidth = Math.max(MIN_BANDWIDTH, 1.06 * std * Math.pow(values.length, -0.2));
  }

  return { xs, density: xs.map((x) => gaussianKdeAt(values, bandwidth, x)) };
}

function quantile(sorted: number[], q: number): number {
  if (sorted.length === 1) return sorted[0];
  const pos = (sorted.length - 1) * q;
  const lo = Math.floor(pos);
  const hi = Math.ceil(pos);
  if (lo === hi) return sorted[lo];
  return sorted[lo] + (sorted[hi] - sorted[lo]) * (pos - lo);
}

// Deterministic pseudo-random jitter (from the point's own id) so dots don't
// reshuffle position on every re-render/refetch.
function jitterFor(id: string): number {
  let hash = 0;
  for (let i = 0; i < id.length; i++) {
    hash = (hash * 31 + id.charCodeAt(i)) | 0;
  }
  return (Math.abs(hash) % 1000) / 1000 - 0.5; // in [-0.5, 0.5]
}

const PLOT_HEIGHT = 260;
const MARGIN = { top: 16, right: 16, bottom: 56, left: 44 };
const SLOT_WIDTH = 150;
const MAX_HALF_WIDTH = 42;
const IQR_BOX_HALF_WIDTH = 8;

function ViolinChart({ metrics, svgRef }: { metrics: CompareMetricResult[]; svgRef: React.RefObject<SVGSVGElement> }) {
  const theme = useTheme();
  const [hovered, setHovered] = React.useState<{ metricId: string; point: CompareDataPoint } | null>(null);

  const plotWidth = SLOT_WIDTH * Math.max(1, metrics.length);
  const width = MARGIN.left + plotWidth + MARGIN.right;
  const height = MARGIN.top + PLOT_HEIGHT + MARGIN.bottom;
  const valueToY = (v: number) => MARGIN.top + (1 - v) * PLOT_HEIGHT;

  // Resolved (non-CSS-var) values -- unlike sx props, colorManipulator's alpha()
  // below can't parse a `var(--...)` reference, so this must read the theme's actual
  // literal colors rather than theme.vars. Same convention WorkspaceDiscovery's
  // similarityColor() already uses for SVG colors, for the same reason.
  const gridColor = theme.palette.divider;
  const axisTextColor = theme.palette.text.secondary;
  const surfaceColor = theme.palette.background.paper;

  return (
    <Box sx={{ overflowX: "auto" }}>
      <svg
        ref={svgRef}
        width={width}
        height={height}
        role="img"
        aria-label="Violin plots comparing selected compounds to the query"
      >
        {/* Gridlines + y-axis labels, shared across every violin (one axis, one scale).
            data-export-role marks colors the SVG export overrides for a white background --
            see chartSvgExport.ts. */}
        {[0, 0.25, 0.5, 0.75, 1].map((tick) => {
          const y = valueToY(tick);
          return (
            <g key={tick}>
              <line
                data-export-role="gridline"
                x1={MARGIN.left}
                x2={width - MARGIN.right}
                y1={y}
                y2={y}
                stroke={gridColor}
                strokeWidth={1}
              />
              <text
                data-export-role="axis-text"
                x={MARGIN.left - 8}
                y={y}
                textAnchor="end"
                dominantBaseline="middle"
                fontSize={11}
                fill={axisTextColor}
              >
                {Math.round(tick * 100)}%
              </text>
            </g>
          );
        })}

        {metrics.map((metric, idx) => {
          const centerX = MARGIN.left + SLOT_WIDTH * (idx + 0.5);
          const values = metric.points.map((p) => p.value).filter((v): v is number => v !== null);
          const sorted = [...values].sort((a, b) => a - b);
          const n = values.length;

          if (n === 0) {
            return (
              <g key={metric.id}>
                <text
                  data-export-role="axis-text"
                  x={centerX}
                  y={MARGIN.top + PLOT_HEIGHT / 2}
                  textAnchor="middle"
                  fontSize={12}
                  fill={axisTextColor}
                >
                  No data
                </text>
                <text
                  data-export-role="axis-text"
                  x={centerX}
                  y={height - MARGIN.bottom + 20}
                  textAnchor="middle"
                  fontSize={12}
                  fill={axisTextColor}
                >
                  {metric.label}
                </text>
              </g>
            );
          }

          const { xs, density } = sampleKde(values);
          const maxDensity = Math.max(...density, 1e-9);
          const scale = MAX_HALF_WIDTH / maxDensity;

          const rightSide = xs.map((x, i) => `${centerX + density[i] * scale},${valueToY(x)}`);
          const leftSide = xs
            .map((x, i) => `${centerX - density[i] * scale},${valueToY(x)}`)
            .reverse();
          const pathD = `M ${rightSide.join(" L ")} L ${leftSide.join(" L ")} Z`;

          const median = quantile(sorted, 0.5);
          const q1 = quantile(sorted, 0.25);
          const q3 = quantile(sorted, 0.75);

          const skipped = metric.points.length - n;

          return (
            <g key={metric.id}>
              <path d={pathD} fill={alpha(metric.color, 0.28)} stroke={metric.color} strokeWidth={2} strokeLinejoin="round" />

              {/* IQR box + median tick */}
              <rect
                x={centerX - IQR_BOX_HALF_WIDTH}
                y={valueToY(q3)}
                width={IQR_BOX_HALF_WIDTH * 2}
                height={Math.max(1, valueToY(q1) - valueToY(q3))}
                fill={alpha(metric.color, 0.35)}
                stroke={metric.color}
                strokeWidth={1}
              />
              <line
                x1={centerX - IQR_BOX_HALF_WIDTH}
                x2={centerX + IQR_BOX_HALF_WIDTH}
                y1={valueToY(median)}
                y2={valueToY(median)}
                stroke={metric.color}
                strokeWidth={2}
              />

              {/* Individual candidates, jittered so overlapping values stay distinguishable */}
              {metric.points.map((point) => {
                if (point.value === null) return null;
                const x = centerX + jitterFor(point.id) * (MAX_HALF_WIDTH * 0.9);
                const y = valueToY(point.value);
                const isHovered = hovered?.metricId === metric.id && hovered.point.id === point.id;
                return (
                  <g
                    key={point.id}
                    onMouseEnter={() => setHovered({ metricId: metric.id, point })}
                    onMouseLeave={() => setHovered((h) => (h?.point.id === point.id ? null : h))}
                    style={{ cursor: "pointer" }}
                  >
                    {/* transparent hit target, bigger than the painted dot */}
                    <circle cx={x} cy={y} r={12} fill="transparent" />
                    <circle
                      data-export-role="dot-ring"
                      cx={x}
                      cy={y}
                      r={isHovered ? 5 : 3.5}
                      fill={metric.color}
                      stroke={surfaceColor}
                      strokeWidth={2}
                    />
                  </g>
                );
              })}

              <text
                data-export-role="axis-text"
                x={centerX}
                y={height - MARGIN.bottom + 20}
                textAnchor="middle"
                fontSize={12}
                fontWeight={500}
                fill={axisTextColor}
              >
                {metric.label}
              </text>
              <text
                data-export-role="axis-text"
                x={centerX}
                y={height - MARGIN.bottom + 36}
                textAnchor="middle"
                fontSize={11}
                fill={axisTextColor}
              >
                {skipped > 0 ? `n=${n} (${skipped} skipped)` : `n=${n}`}
              </text>
            </g>
          );
        })}
      </svg>

      {hovered && (
        <Typography variant="caption" color="text.secondary" sx={{ display: "block", mt: 0.5 }}>
          <strong>{hovered.point.label}</strong> — {METRIC_LABELS[hovered.metricId as CompareMetricId] ?? hovered.metricId}:{" "}
          {((hovered.point.value ?? 0) * 100).toFixed(1)}%
        </Typography>
      )}
    </Box>
  );
}

// ---- Panel: metric picker + data wiring ------------------------------------

function chartColorFor(theme: Theme, metricId: CompareMetricId): string {
  // Resolved values, not theme.vars -- these feed alpha() in ViolinChart, which
  // can't parse a `var(--...)` CSS variable reference.
  const { palette } = theme;
  const colors = [palette.primary.main, palette.warning.main, palette.success.main, palette.error.main];
  // Slot fixed by position in METRIC_DEFS (not by which metrics happen to be
  // enabled), so a metric's color stays its identity every time it's toggled.
  const slot = METRIC_DEFS.findIndex((def) => def.id === metricId);
  return colors[Math.max(0, slot) % colors.length];
}

export function WorkspaceCompare({
  results,
  queryOriginSmiles,
}: {
  results: DiscoveryResult[];
  queryOriginSmiles: string | null;
}) {
  const theme = useTheme();

  const [enabledMetrics, setEnabledMetrics] = React.useState<Set<CompareMetricId>>(
    new Set<CompareMetricId>(["fingerprintCosine", "normalizedAlignment"])
  );
  const [radius, setRadius] = React.useState<number>(MORGAN_RADIUS_DEFAULT);
  const [nBits, setNBits] = React.useState<number>(MORGAN_NBITS_DEFAULT);
  const [pickerAnchor, setPickerAnchor] = React.useState<HTMLElement | null>(null);
  const svgRef = React.useRef<SVGSVGElement>(null);

  const toggleMetric = (id: CompareMetricId) => {
    setEnabledMetrics((prev) => {
      const next = new Set(prev);
      if (next.has(id)) next.delete(id);
      else next.add(id);
      return next;
    });
  };

  const tanimotoEnabled = enabledMetrics.has("tanimoto");
  const tanimotoCandidates = React.useMemo(
    () =>
      results
        .filter((r): r is DiscoveryResult & { smiles: string } => r.type === "compound" && !!r.smiles)
        .map((r) => ({ id: r.entryId, smiles: r.smiles })),
    [results]
  );

  const tanimotoQuery = useQuery({
    queryKey: [
      "discoveryCompareTanimoto",
      queryOriginSmiles,
      radius,
      nBits,
      tanimotoCandidates.map((c) => c.id).join(","),
    ],
    queryFn: ({ signal }) =>
      runDiscoveryCompareTanimoto(
        { querySmiles: queryOriginSmiles as string, radius, nBits, entries: tanimotoCandidates },
        signal
      ),
    enabled: tanimotoEnabled && !!queryOriginSmiles && tanimotoCandidates.length > 0,
    // Result is a pure function of the query key (query smiles + candidate smiles +
    // radius/nBits) -- never go stale, so toggling the metric off/on or switching
    // Discovery views and back doesn't re-run the RDKit fingerprinting server-side.
    staleTime: Infinity,
  });

  const metrics: CompareMetricResult[] = [];

  if (enabledMetrics.has("fingerprintCosine")) {
    metrics.push({
      id: "fingerprintCosine",
      label: METRIC_LABELS.fingerprintCosine,
      color: chartColorFor(theme, "fingerprintCosine"),
      points: results.map((r) => ({ id: r.entryId, label: r.name, value: r.fingerprintSimilarity })),
    });
  }

  if (enabledMetrics.has("normalizedAlignment")) {
    metrics.push({
      id: "normalizedAlignment",
      label: METRIC_LABELS.normalizedAlignment,
      color: chartColorFor(theme, "normalizedAlignment"),
      points: results.map((r) => ({ id: r.entryId, label: r.name, value: r.normalizedAlignmentScorePct / 100 })),
    });
  }

  if (tanimotoEnabled && queryOriginSmiles && tanimotoCandidates.length > 0) {
    const valueById = new Map((tanimotoQuery.data?.values ?? []).map((v) => [v.id, v.tanimotoSimilarity]));
    const nameById = new Map(results.map((r) => [r.entryId, r.name]));
    metrics.push({
      id: "tanimoto",
      label: `Tanimoto (r${radius}, ${nBits}b)`,
      color: chartColorFor(theme, "tanimoto"),
      points: tanimotoCandidates.map((c) => ({
        id: c.id,
        label: nameById.get(c.id) ?? c.id,
        value: valueById.get(c.id) ?? null,
      })),
      isLoading: tanimotoQuery.isFetching,
      error: tanimotoQuery.error ? (tanimotoQuery.error as Error).message : null,
    });
  }

  return (
    <Box>
      <Typography variant="caption" color="text.secondary" sx={{ display: "block", mb: 1.5 }}>
        Each violin shows the distribution of one metric across the selected results, comparing every candidate back
        to the query. Pick which metrics to include below.
      </Typography>

      {/* A button + Popover rather than one control per metric -- the toolbar stays
          exactly this wide no matter how many metrics METRIC_DEFS grows to. The
          enabled set is still visible at a glance via the chips beside the button. */}
      <Stack direction="row" spacing={1} alignItems="center" flexWrap="wrap" sx={{ mb: 1.5 }}>
        <Button
          size="small"
          variant="outlined"
          startIcon={<TuneIcon fontSize="small" />}
          endIcon={<ArrowDropDownIcon />}
          onClick={(e) => setPickerAnchor(e.currentTarget)}
        >
          Metrics{enabledMetrics.size > 0 ? ` (${enabledMetrics.size})` : ""}
        </Button>

        {enabledMetrics.size === 0 ? (
          <Typography variant="caption" color="text.secondary">
            None selected
          </Typography>
        ) : (
          METRIC_DEFS.filter((def) => enabledMetrics.has(def.id)).map((def) => (
            <Chip
              key={def.id}
              size="small"
              variant="outlined"
              label={def.id === "tanimoto" ? `${def.label} (r${radius}, ${nBits}b)` : def.label}
              onDelete={() => toggleMetric(def.id)}
              sx={{
                borderColor: chartColorFor(theme, def.id),
                color: chartColorFor(theme, def.id),
              }}
            />
          ))
        )}
      </Stack>

      <Popover
        open={Boolean(pickerAnchor)}
        anchorEl={pickerAnchor}
        onClose={() => setPickerAnchor(null)}
        anchorOrigin={{ vertical: "bottom", horizontal: "left" }}
        transformOrigin={{ vertical: "top", horizontal: "left" }}
      >
        <Box sx={{ p: 2, width: 320 }}>
          <Typography variant="subtitle2" sx={{ mb: 1 }}>
            Metrics to plot
          </Typography>
          <Stack divider={<Divider flexItem />} spacing={1}>
            {METRIC_DEFS.map((def) => (
              <Box key={def.id} sx={{ pb: 0.5 }}>
                <FormControlLabel
                  sx={{ alignItems: "flex-start", m: 0 }}
                  control={
                    <Checkbox
                      size="small"
                      sx={{ pt: 0 }}
                      checked={enabledMetrics.has(def.id)}
                      onChange={() => toggleMetric(def.id)}
                    />
                  }
                  label={
                    <Box>
                      <Typography variant="body2">{def.label}</Typography>
                      <Typography variant="caption" color="text.secondary">
                        {def.description}
                      </Typography>
                    </Box>
                  }
                />

                {def.id === "tanimoto" && tanimotoEnabled && (
                  <Stack direction="row" spacing={1} sx={{ pl: 4.5, mt: 1 }}>
                    <TextField
                      select
                      size="small"
                      label="Radius"
                      value={radius}
                      onChange={(e) => setRadius(Number(e.target.value))}
                      sx={{ width: 100 }}
                    >
                      {[1, 2, 3, 4].map((r) => (
                        <MenuItem key={r} value={r}>
                          {r}
                        </MenuItem>
                      ))}
                    </TextField>
                    <TextField
                      select
                      size="small"
                      label="Fingerprint size"
                      value={nBits}
                      onChange={(e) => setNBits(Number(e.target.value))}
                      sx={{ width: 140 }}
                    >
                      {MORGAN_NBITS_OPTIONS.map((b) => (
                        <MenuItem key={b} value={b}>
                          {b} bits
                        </MenuItem>
                      ))}
                    </TextField>
                  </Stack>
                )}
              </Box>
            ))}
          </Stack>
        </Box>
      </Popover>

      {tanimotoEnabled && !queryOriginSmiles && (
        <Alert severity="info" sx={{ mb: 1.5 }}>
          Tanimoto similarity compares real molecular structures, so it needs the query's own SMILES. Pick an
          uploaded compound above and click "Use this" before querying to enable it.
        </Alert>
      )}
      {tanimotoEnabled && queryOriginSmiles && tanimotoCandidates.length === 0 && (
        <Alert severity="info" sx={{ mb: 1.5 }}>
          None of the selected results are compounds with a known structure (SMILES), so there's nothing to compare
          for this metric.
        </Alert>
      )}
      {tanimotoEnabled && tanimotoQuery.error && (
        <Alert severity="error" sx={{ mb: 1.5 }}>
          {(tanimotoQuery.error as Error).message || "Failed to compute Tanimoto similarity."}
        </Alert>
      )}

      {metrics.length === 0 ? (
        <Typography variant="body2" color="text.secondary">
          Pick at least one metric above to plot.
        </Typography>
      ) : (
        <Box>
          {metrics.some((m) => m.isLoading) && (
            <Stack direction="row" spacing={1} alignItems="center" sx={{ mb: 1 }}>
              <CircularProgress size={14} />
              <Typography variant="caption" color="text.secondary">
                Recalculating…
              </Typography>
            </Stack>
          )}
          <ViolinChart metrics={metrics} svgRef={svgRef} />
          <Button
            size="small"
            variant="text"
            startIcon={<DownloadIcon fontSize="small" />}
            onClick={() => svgRef.current && exportChartSvg(svgRef.current, "compound_comparison.svg")}
            sx={{ mt: 0.5 }}
          >
            Download SVG
          </Button>
        </Box>
      )}
    </Box>
  );
}
