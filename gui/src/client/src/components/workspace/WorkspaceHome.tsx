import React from "react";
import Box from "@mui/material/Box";
import Card from "@mui/material/Card";
import CardContent from "@mui/material/CardContent";
import MuiLink from "@mui/material/Link";
import Typography from "@mui/material/Typography";
import Alert from "@mui/material/Alert";
import Skeleton from "@mui/material/Skeleton";
import { useTheme } from "@mui/material/styles";
import { Link as RouterLink } from "react-router-dom";
import { PieChart } from "@mui/x-charts/PieChart";
import { BarChart } from "@mui/x-charts/BarChart";
import { getDatabaseStats, getAnnotationStats } from "../../features/database/api";
import { DatabaseStatsResp, AnnotationStatsResp, Count } from "../../features/database/types";

const ENTRY_TYPE_LABELS: Record<string, string> = {
  compound: "Compounds",
  bgc: "Gene clusters (BGCs)",
};

const StatTile: React.FC<{ label: string; value: string; caption?: string }> = ({ label, value, caption }) => (
  <Card variant="outlined" sx={{ flex: "1 1 180px", minWidth: 180 }}>
    <CardContent>
      <Typography variant="body2" color="text.secondary">
        {label}
      </Typography>
      <Typography variant="h4" component="div" sx={{ fontWeight: 600, lineHeight: 1.3 }}>
        {value}
      </Typography>
      {caption && (
        <Typography variant="caption" color="text.secondary">
          {caption}
        </Typography>
      )}
    </CardContent>
  </Card>
);

const ChartCard: React.FC<{ title: string; description?: string; children: React.ReactNode }> = ({
  title,
  description,
  children,
}) => (
  <Card variant="outlined" sx={{ flex: "1 1 380px", minWidth: 0 }}>
    <CardContent>
      <Typography component="h2" variant="subtitle1" sx={{ mb: description ? 0.25 : 1 }}>
        {title}
      </Typography>
      {description && (
        <Typography variant="caption" color="text.secondary" sx={{ display: "block", mb: 1 }}>
          {description}
        </Typography>
      )}
      {children}
    </CardContent>
  </Card>
);

const AnnotationBarCard: React.FC<{
  title: string;
  description?: string;
  counts: Count[];
  color: string;
  emptyMessage: string;
}> = ({ title, description, counts, color, emptyMessage }) => (
  <ChartCard title={title} description={description}>
    {counts.length === 0 ? (
      <Typography variant="body2" color="text.secondary">
        {emptyMessage}
      </Typography>
    ) : (
      <BarChart
        dataset={counts.map((c) => ({ label: c.label, count: c.count }))}
        xAxis={[{ scaleType: "band", dataKey: "label" }]}
        series={[{ dataKey: "count", color }]}
        height={220}
      />
    )}
  </ChartCard>
);

export const WorkspaceHome: React.FC = () => {
  const theme = useTheme();
  const palette = theme.vars || theme;

  const [stats, setStats] = React.useState<DatabaseStatsResp | null>(null);
  const [loading, setLoading] = React.useState<boolean>(true);
  const [error, setError] = React.useState<string | null>(null);

  const [annotationStats, setAnnotationStats] = React.useState<AnnotationStatsResp | null>(null);
  const [annotationLoading, setAnnotationLoading] = React.useState<boolean>(true);
  const [annotationError, setAnnotationError] = React.useState<string | null>(null);

  React.useEffect(() => {
    const controller = new AbortController();
    setLoading(true);
    setError(null);

    getDatabaseStats(controller.signal)
      .then(setStats)
      .catch((err) => {
        if (controller.signal.aborted) return;
        setError(err instanceof Error ? err.message : String(err));
      })
      .finally(() => {
        if (controller.signal.aborted) return;
        setLoading(false);
      });

    return () => controller.abort();
  }, []);

  React.useEffect(() => {
    const controller = new AbortController();
    setAnnotationLoading(true);
    setAnnotationError(null);

    getAnnotationStats(controller.signal)
      .then(setAnnotationStats)
      .catch((err) => {
        if (controller.signal.aborted) return;
        setAnnotationError(err instanceof Error ? err.message : String(err));
      })
      .finally(() => {
        if (controller.signal.aborted) return;
        setAnnotationLoading(false);
      });

    return () => controller.abort();
  }, []);

  const chartColors = [
    palette.palette.primary.main,
    palette.palette.warning.main,
    palette.palette.success.main,
    palette.palette.error.main,
    palette.palette.info.main,
    palette.palette.secondary.main,
  ];

  return (
    <Box
      sx={{
        width: "100%",
        mx: "auto",
        display: "flex",
        flexDirection: "column",
        gap: "16px",
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
            Use the navigation menu on the left to upload your data in the&nbsp;
            <MuiLink
              component={RouterLink}
              to="/dashboard/upload"
              underline="hover"
              color={(theme.vars || theme).palette.primary.main}
              sx={{ fontWeight: "500" }}
            >
              Upload tab
            </MuiLink>
            .
            Biosynthetic fingerprints are automatically parsed from your data during upload.
            You can use the biosynthetic fingerprints parsed from you data for exploratory data analysis in the&nbsp;
            <MuiLink
              component={RouterLink}
              to="/dashboard/discovery"
              underline="hover"
              color={(theme.vars || theme).palette.primary.main}
              sx={{ fontWeight: "500" }}
            >
              Discovery tab
            </MuiLink>
            .
            Exploratory data analysis allows you to retrieve similar biosynthetic fingerprints and their associated metadata from the database using your uploaded data.
          </Typography>
        </CardContent>
      </Card>

      <Typography component="h2" variant="subtitle1" sx={{ mt: 1 }}>
        Database at a glance
      </Typography>

      {error && (
        <Alert severity="warning" variant="outlined">
          Couldn't load database statistics: {error}
        </Alert>
      )}

      {!error && loading && (
        <Box sx={{ display: "flex", flexWrap: "wrap", gap: "16px" }}>
          {Array.from({ length: 3 }).map((_, i) => (
            <Skeleton key={i} variant="rounded" height={98} sx={{ flex: "1 1 180px", minWidth: 180 }} />
          ))}
        </Box>
      )}

      {!error && !loading && stats && (
        <>
          <Box sx={{ display: "flex", flexWrap: "wrap", gap: "16px" }}>
            <StatTile label="Entries in the database" value={stats.totalEntries.toLocaleString()} />
            <StatTile
              label="Unique building blocks"
              value={stats.uniqueBlockCount.toLocaleString()}
              caption="Monomers & tailoring events"
            />
            <StatTile
              label="Avg. sequence length"
              value={stats.sequenceLengthAvg.toFixed(1)}
              caption={`Range ${stats.sequenceLengthMin}–${stats.sequenceLengthMax} blocks`}
            />
          </Box>

          <Box sx={{ display: "flex", flexWrap: "wrap", gap: "16px" }}>
            <ChartCard title="Entry types" description="What kind of data each entry represents">
              <PieChart
                series={[
                  {
                    data: stats.countsByType.map((c, i) => ({
                      id: i,
                      value: c.count,
                      label: ENTRY_TYPE_LABELS[c.label] ?? c.label,
                    })),
                    innerRadius: 45,
                    paddingAngle: 2,
                    cornerRadius: 3,
                    highlightScope: { fade: "global", highlight: "item" },
                  },
                ]}
                colors={chartColors}
                height={220}
                slotProps={{ legend: { direction: "vertical" } }}
              />
            </ChartCard>
          </Box>

          <Typography component="h2" variant="subtitle1" sx={{ mt: 1 }}>
            Annotations
          </Typography>

          {annotationError && (
            <Alert severity="warning" variant="outlined">
              Couldn't load annotation statistics: {annotationError}
            </Alert>
          )}

          {!annotationError && annotationLoading && (
            <Box sx={{ display: "flex", flexWrap: "wrap", gap: "16px" }}>
              {Array.from({ length: 2 }).map((_, i) => (
                <Skeleton key={i} variant="rounded" height={260} sx={{ flex: "1 1 380px", minWidth: 0 }} />
              ))}
            </Box>
          )}

          {!annotationError && !annotationLoading && annotationStats && (
            <>
              <Box sx={{ display: "flex", flexWrap: "wrap", gap: "16px" }}>
                <StatTile
                  label="Entries with annotations"
                  value={annotationStats.withAnnotationCount.toLocaleString()}
                  caption={`${annotationStats.withoutAnnotationCount.toLocaleString()} without`}
                />
              </Box>

              <Typography component="h3" variant="subtitle2" sx={{ mt: 1 }}>
                Phylogeny
              </Typography>
              <Box sx={{ display: "flex", flexWrap: "wrap", gap: "16px" }}>
                <AnnotationBarCard
                  title="Type"
                  description="Bacterium / fungus / archaeon / other"
                  counts={annotationStats.phylogenyTypeCounts}
                  color={chartColors[0]}
                  emptyMessage="No phylogeny annotations yet."
                />
                <AnnotationBarCard
                  title="Genus"
                  description="Most common genera"
                  counts={annotationStats.phylogenyGenusCounts}
                  color={chartColors[1]}
                  emptyMessage="No phylogeny annotations yet."
                />
                <AnnotationBarCard
                  title="Species"
                  description="Most common species"
                  counts={annotationStats.phylogenySpeciesCounts}
                  color={chartColors[2]}
                  emptyMessage="No phylogeny annotations yet."
                />
              </Box>

              <Typography component="h3" variant="subtitle2" sx={{ mt: 1 }}>
                Chemical class
              </Typography>
              <Typography variant="caption" color="text.secondary" sx={{ display: "block", mt: -1 }}>
                NPClassifier, predicted from every compound's own structure
              </Typography>
              <Box sx={{ display: "flex", flexWrap: "wrap", gap: "16px" }}>
                <AnnotationBarCard
                  title="Pathway"
                  counts={annotationStats.chemicalClassPathwayCounts}
                  color={chartColors[3]}
                  emptyMessage="No NPClassifier annotations yet."
                />
                <AnnotationBarCard
                  title="Superclass"
                  counts={annotationStats.chemicalClassSuperclassCounts}
                  color={chartColors[4]}
                  emptyMessage="No NPClassifier annotations yet."
                />
                <AnnotationBarCard
                  title="Class"
                  counts={annotationStats.chemicalClassClassCounts}
                  color={chartColors[5]}
                  emptyMessage="No NPClassifier annotations yet."
                />
              </Box>

              <Typography component="h3" variant="subtitle2" sx={{ mt: 1 }}>
                Biosynthetic class
              </Typography>
              <Typography variant="caption" color="text.secondary" sx={{ display: "block", mt: -1 }}>
                MIBiG's own coarse label (PKS / NRPS / RiPP / ...), a separate classification from NPClassifier's chemical class above
              </Typography>
              <Box sx={{ display: "flex", flexWrap: "wrap", gap: "16px" }}>
                <AnnotationBarCard
                  title="Biosynthetic class"
                  counts={annotationStats.biosyntheticClassCounts}
                  color={chartColors[0]}
                  emptyMessage="No biosynthetic class annotations yet."
                />
              </Box>

              <Typography component="h3" variant="subtitle2" sx={{ mt: 1 }}>
                Bioactivity
              </Typography>
              <Typography variant="caption" color="text.secondary" sx={{ display: "block", mt: -1 }}>
                ChEMBL + ChEBI, looked up by structure (InChIKey)
              </Typography>
              <Box sx={{ display: "flex", flexWrap: "wrap", gap: "16px" }}>
                <AnnotationBarCard
                  title="ATC category"
                  description="WHO therapeutic classification (ChEMBL)"
                  counts={annotationStats.bioactivityAtcCounts}
                  color={chartColors[1]}
                  emptyMessage="No ChEMBL ATC matches yet."
                />
                <AnnotationBarCard
                  title="Clinical phase"
                  description="Furthest development stage reached (ChEMBL)"
                  counts={annotationStats.bioactivityMaxPhaseCounts}
                  color={chartColors[2]}
                  emptyMessage="No ChEMBL clinical-phase matches yet."
                />
                <AnnotationBarCard
                  title="Biological role"
                  description="ChEBI role ontology"
                  counts={annotationStats.bioactivityBiologicalRoleCounts}
                  color={chartColors[3]}
                  emptyMessage="No ChEBI biological-role matches yet."
                />
                <AnnotationBarCard
                  title="Chemical role"
                  description="ChEBI role ontology"
                  counts={annotationStats.bioactivityChemicalRoleCounts}
                  color={chartColors[4]}
                  emptyMessage="No ChEBI chemical-role matches yet."
                />
              </Box>
            </>
          )}
        </>
      )}
    </Box>
  );
};
