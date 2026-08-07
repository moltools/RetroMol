import React from "react";
import Alert from "@mui/material/Alert";
import Box from "@mui/material/Box";
import Button from "@mui/material/Button";
import CircularProgress from "@mui/material/CircularProgress";
import Stack from "@mui/material/Stack";
import Typography from "@mui/material/Typography";
import DownloadIcon from "@mui/icons-material/Download";
import { useQuery } from "@tanstack/react-query";
import { useTheme } from "@mui/material/styles";
import { runDiscoveryShape } from "../../features/discovery/api";
import { MAX_SHAPE_ENTRIES, type DiscoveryResult } from "../../features/discovery/types";
import { exportChartSvg } from "./chartSvgExport";

type ShapePointKind = "query" | "candidate";
type ShapeOrigin = "database" | "upload";

type ShapePoint = {
  id: string;
  label: string;
  npr1: number;
  npr2: number;
  kind: ShapePointKind;
  origin: ShapeOrigin;
};

// Triangle envelope corners in (NPR1, NPR2) space -- identical to the reference
// script's triangle_x/triangle_y: rod (0,1), disc (0.5,0.5), sphere (1,1).
const CORNER_ROD: [number, number] = [0, 1];
const CORNER_DISC: [number, number] = [0.5, 0.5];
const CORNER_SPHERE: [number, number] = [1, 1];

// Padded a bit beyond the triangle itself so corner labels have room; kept at equal
// px-per-unit on both axes (see plotHeight below) so the triangle isn't stretched.
const X_DOMAIN: [number, number] = [-0.1, 1.1];
const Y_DOMAIN: [number, number] = [0.4, 1.1];
const PLOT_WIDTH = 380;
const MARGIN = { top: 22, right: 22, bottom: 26, left: 22 };
const LEGEND_HEIGHT = 28;
const POINT_RADIUS = 4.5;
const QUERY_MARKER_HALF_SIZE = 6.5;

function ShapeChart({ points, svgRef }: { points: ShapePoint[]; svgRef: React.RefObject<SVGSVGElement> }) {
  const theme = useTheme();
  const [hovered, setHovered] = React.useState<ShapePoint | null>(null);

  const pxPerUnit = PLOT_WIDTH / (X_DOMAIN[1] - X_DOMAIN[0]);
  const plotHeight = pxPerUnit * (Y_DOMAIN[1] - Y_DOMAIN[0]);
  const width = MARGIN.left + PLOT_WIDTH + MARGIN.right;
  const height = MARGIN.top + plotHeight + MARGIN.bottom + LEGEND_HEIGHT;

  const xToPx = (x: number) => MARGIN.left + (x - X_DOMAIN[0]) * pxPerUnit;
  const yToPx = (y: number) => MARGIN.top + (Y_DOMAIN[1] - y) * pxPerUnit;

  // Resolved (non-CSS-var) values -- these feed straight into SVG fill/stroke
  // attributes, which can't resolve `var(--...)` the way sx props can. Same
  // convention as ViolinChart/WorkspaceCompare.
  const envelopeColor = theme.palette.text.secondary;
  const axisTextColor = theme.palette.text.secondary;
  const surfaceColor = theme.palette.background.paper;
  // primary/warning/success -- same categorical rotation WorkspaceCompare's
  // chartColorFor uses, chosen so database/query/upload read as three clearly
  // distinct hues (info.main is too close to primary.main's blue to tell apart).
  const databaseColor = theme.palette.primary.main;
  const uploadColor = theme.palette.success.main;
  const queryColor = theme.palette.warning.main;

  const pointColor = (p: ShapePoint) => (p.kind === "query" ? queryColor : p.origin === "upload" ? uploadColor : databaseColor);

  const hasDatabase = points.some((p) => p.kind === "candidate" && p.origin === "database");
  const hasUpload = points.some((p) => p.kind === "candidate" && p.origin === "upload");
  const hasQuery = points.some((p) => p.kind === "query");

  const legendItems: { label: string; color: string; shape: "circle" | "diamond" }[] = [
    ...(hasQuery ? [{ label: "Query", color: queryColor, shape: "diamond" as const }] : []),
    ...(hasDatabase ? [{ label: "Database", color: databaseColor, shape: "circle" as const }] : []),
    ...(hasUpload ? [{ label: "Upload", color: uploadColor, shape: "circle" as const }] : []),
  ];

  const trianglePath = `M ${xToPx(CORNER_ROD[0])},${yToPx(CORNER_ROD[1])} L ${xToPx(CORNER_DISC[0])},${yToPx(
    CORNER_DISC[1]
  )} L ${xToPx(CORNER_SPHERE[0])},${yToPx(CORNER_SPHERE[1])} Z`;

  const legendY = MARGIN.top + plotHeight + MARGIN.bottom + LEGEND_HEIGHT / 2;

  return (
    <Box sx={{ overflowX: "auto" }}>
      <svg ref={svgRef} width={width} height={height} role="img" aria-label="PMI shape triangle: rod, disc, and sphere character of each compound">
        <path
          data-export-role="gridline"
          d={trianglePath}
          fill="none"
          stroke={envelopeColor}
          strokeWidth={1.5}
          strokeDasharray="5,4"
        />

        <text data-export-role="axis-text" x={xToPx(0)} y={yToPx(1) - 8} textAnchor="middle" fontSize={12} fontWeight={600} fill={axisTextColor}>
          Rod
        </text>
        <text data-export-role="axis-text" x={xToPx(1)} y={yToPx(1) - 8} textAnchor="middle" fontSize={12} fontWeight={600} fill={axisTextColor}>
          Sphere
        </text>
        <text data-export-role="axis-text" x={xToPx(0.5)} y={yToPx(0.5) + 18} textAnchor="middle" fontSize={12} fontWeight={600} fill={axisTextColor}>
          Disc
        </text>

        <text data-export-role="axis-text" x={xToPx(0.5)} y={height - 6} textAnchor="middle" fontSize={11} fill={axisTextColor}>
          NPR1 (I₁/I₃) →
        </text>
        <text
          data-export-role="axis-text"
          x={12}
          y={MARGIN.top + plotHeight / 2}
          textAnchor="middle"
          fontSize={11}
          fill={axisTextColor}
          transform={`rotate(-90, 12, ${MARGIN.top + plotHeight / 2})`}
        >
          NPR2 (I₂/I₃) →
        </text>

        {points.map((p) => {
          const cx = xToPx(p.npr1);
          const cy = yToPx(p.npr2);
          const color = pointColor(p);
          const isHovered = hovered?.id === p.id;
          const r = isHovered ? POINT_RADIUS + 1.5 : POINT_RADIUS;

          return (
            <g
              key={p.id}
              onMouseEnter={() => setHovered(p)}
              onMouseLeave={() => setHovered((h) => (h?.id === p.id ? null : h))}
              style={{ cursor: "pointer" }}
            >
              {/* transparent hit target, bigger than the painted mark */}
              <circle cx={cx} cy={cy} r={14} fill="transparent" />
              {p.kind === "query" ? (
                <polygon
                  data-export-role="dot-ring"
                  points={`${cx},${cy - QUERY_MARKER_HALF_SIZE} ${cx + QUERY_MARKER_HALF_SIZE},${cy} ${cx},${
                    cy + QUERY_MARKER_HALF_SIZE
                  } ${cx - QUERY_MARKER_HALF_SIZE},${cy}`}
                  fill={color}
                  stroke={surfaceColor}
                  strokeWidth={2}
                  strokeLinejoin="round"
                />
              ) : (
                <circle data-export-role="dot-ring" cx={cx} cy={cy} r={r} fill={color} stroke={surfaceColor} strokeWidth={2} />
              )}
            </g>
          );
        })}

        {legendItems.map((item, idx) => {
          const itemX = MARGIN.left + idx * 100;
          return (
            <g key={item.label}>
              {item.shape === "diamond" ? (
                <polygon
                  points={`${itemX + 6},${legendY - 6} ${itemX + 12},${legendY} ${itemX + 6},${legendY + 6} ${itemX},${legendY}`}
                  fill={item.color}
                />
              ) : (
                <circle cx={itemX + 6} cy={legendY} r={5} fill={item.color} />
              )}
              <text data-export-role="axis-text" x={itemX + 18} y={legendY} dominantBaseline="middle" fontSize={11} fill={axisTextColor}>
                {item.label}
              </text>
            </g>
          );
        })}
      </svg>

      {hovered && (
        <Typography variant="caption" color="text.secondary" sx={{ display: "block", mt: 0.5 }}>
          <strong>{hovered.label}</strong> — NPR1: {hovered.npr1.toFixed(3)}, NPR2: {hovered.npr2.toFixed(3)}
        </Typography>
      )}
    </Box>
  );
}

type PrecomputedShape = {
  results: { id: string; npr1: number; npr2: number }[];
  skipped: { id: string; reason: string }[];
};

export function WorkspaceShape({
  results,
  queryOriginSmiles,
  precomputedShape = null,
}: {
  results: DiscoveryResult[];
  queryOriginSmiles: string | null;
  // Precomputed at query-submit time (see DialogViewDiscoveryQuery). When present,
  // used directly instead of hitting /api/discoveryShape -- this view has no
  // interactive conformer-budget picker, so unlike WorkspaceCompare's Tanimoto
  // metric there's no "still matches" check needed, it's simply always preferred.
  precomputedShape?: PrecomputedShape | null;
}) {
  const svgRef = React.useRef<SVGSVGElement>(null);

  const shapeCandidates = React.useMemo(
    () =>
      results
        .filter((r): r is DiscoveryResult & { smiles: string } => r.type === "compound" && !!r.smiles)
        .map((r) => ({ id: r.entryId, smiles: r.smiles, name: r.name, origin: r.origin })),
    [results]
  );

  const totalEntries = shapeCandidates.length + (queryOriginSmiles ? 1 : 0);
  const tooManyEntries = totalEntries > MAX_SHAPE_ENTRIES;

  const shapeQuery = useQuery({
    queryKey: ["discoveryShape", queryOriginSmiles, shapeCandidates.map((c) => c.id).join(",")],
    queryFn: ({ signal }) =>
      runDiscoveryShape(
        {
          entries: [
            ...(queryOriginSmiles ? [{ id: "query", smiles: queryOriginSmiles }] : []),
            ...shapeCandidates.map((c) => ({ id: c.id, smiles: c.smiles })),
          ],
        },
        signal
      ),
    enabled: shapeCandidates.length > 0 && !tooManyEntries && !precomputedShape,
    // This computation is expensive (conformer search + MMFF optimization per
    // molecule) and its result is a pure function of the query key (selection +
    // conformer budget) -- never go stale, so flipping to another Discovery view
    // and back doesn't silently re-run it against the server. A genuinely new
    // computation only happens when the key itself changes (selection edited).
    staleTime: Infinity,
  });

  const nameById = React.useMemo(() => new Map(shapeCandidates.map((c) => [c.id, c.name])), [shapeCandidates]);
  const originById = React.useMemo(() => new Map(shapeCandidates.map((c) => [c.id, c.origin])), [shapeCandidates]);

  // Precomputed data was captured over every top-X result at query-submit time, not
  // just the currently selected subset -- filter down to the current selection (plus
  // the query itself) so it behaves like the live-fetch path, which only ever
  // requested exactly `shapeCandidates` in the first place.
  const candidateIds = React.useMemo(() => new Set(shapeCandidates.map((c) => c.id)), [shapeCandidates]);
  const inSelection = (id: string) => id === "query" || candidateIds.has(id);

  const points: ShapePoint[] = (
    precomputedShape ? precomputedShape.results.filter((r) => inSelection(r.id)) : shapeQuery.data?.results ?? []
  ).map((r) => ({
    id: r.id,
    label: r.id === "query" ? "Query" : nameById.get(r.id) ?? r.id,
    npr1: r.npr1,
    npr2: r.npr2,
    kind: r.id === "query" ? "query" : "candidate",
    origin: (originById.get(r.id) ?? "database") as ShapeOrigin,
  }));

  const skipped = precomputedShape
    ? precomputedShape.skipped.filter((s) => inSelection(s.id))
    : shapeQuery.data?.skipped ?? [];

  return (
    <Box>
      <Typography variant="caption" color="text.secondary" sx={{ display: "block", mb: 1 }}>
        Classifies each compound's overall 3D shape as rod-, disc-, or sphere-like from its lowest-energy conformer's
        normalized principal moments of inertia (Sauer &amp; Schwarz, <em>J. Chem. Inf. Comput. Sci.</em> 2003).
      </Typography>

      <Alert severity="info" sx={{ mb: 1.5 }}>
        This runs a reduced-cost conformer search (~{25} conformers) so it stays responsive in the browser, and not the
        100+ conformers a literature-quality analysis typically uses. Treat these shapes as approximate; for
        publication-quality PMI values, run the app locally with a larger conformer budget.
      </Alert>

      {shapeCandidates.length === 0 ? (
        <Typography variant="body2" color="text.secondary">
          None of the selected results are compounds with a known structure (SMILES), so there's nothing to plot.
        </Typography>
      ) : tooManyEntries ? (
        <Alert severity="warning" sx={{ mb: 1.5 }}>
          {totalEntries} compounds selected -- shape analysis is limited to {MAX_SHAPE_ENTRIES} at a time. Narrow
          your selection above to compute shapes.
        </Alert>
      ) : (
        <>
          {shapeQuery.isLoading ? (
            <CircularProgress size={20} />
          ) : shapeQuery.error ? (
            <Alert severity="error">{(shapeQuery.error as Error).message || "Failed to compute molecular shapes."}</Alert>
          ) : (
            <Box>
              {shapeQuery.isFetching && (
                <Stack direction="row" spacing={1} alignItems="center" sx={{ mb: 1 }}>
                  <CircularProgress size={14} />
                  <Typography variant="caption" color="text.secondary">
                    Recalculating…
                  </Typography>
                </Stack>
              )}
              {skipped.length > 0 && (
                <Alert severity="warning" sx={{ mb: 1.5 }}>
                  Skipped {skipped.length} compound{skipped.length === 1 ? "" : "s"}: {skipped[0].reason}
                  {skipped.length > 1 ? ` (and ${skipped.length - 1} more)` : ""}
                </Alert>
              )}
              <ShapeChart points={points} svgRef={svgRef} />
              <Button
                size="small"
                variant="text"
                startIcon={<DownloadIcon fontSize="small" />}
                onClick={() => svgRef.current && exportChartSvg(svgRef.current, "pmi_shape_triangle.svg")}
                sx={{ mt: 0.5 }}
              >
                Download SVG
              </Button>
            </Box>
          )}
        </>
      )}
    </Box>
  );
}
