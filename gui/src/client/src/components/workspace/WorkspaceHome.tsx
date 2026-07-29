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
import { getDatabaseStats } from "../../features/database/api";
import { DatabaseStatsResp } from "../../features/database/types";

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

export const WorkspaceHome: React.FC = () => {
  const theme = useTheme();
  const palette = theme.vars || theme;

  const [stats, setStats] = React.useState<DatabaseStatsResp | null>(null);
  const [loading, setLoading] = React.useState<boolean>(true);
  const [error, setError] = React.useState<string | null>(null);

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

  const chartColors = [
    palette.palette.primary.main,
    palette.palette.warning.main,
    palette.palette.success.main,
    palette.palette.error.main,
    palette.palette.info.main,
    palette.palette.secondary.main,
  ];

  const resolvedPct = stats && stats.totalEntries > 0 ? (100 * stats.fullyResolvedCount) / stats.totalEntries : 0;

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
          {Array.from({ length: 4 }).map((_, i) => (
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
            <StatTile
              label="Fully resolved sequences"
              value={`${resolvedPct.toFixed(0)}%`}
              caption="No unidentified blocks"
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

            <ChartCard title="Sequence resolution" description="Share of entries RetroMol could fully identify">
              <PieChart
                series={[
                  {
                    data: [
                      { id: 0, value: stats.fullyResolvedCount, label: "Fully resolved" },
                      { id: 1, value: stats.hasUnknownBlockCount, label: "Contains unresolved block" },
                    ],
                    innerRadius: 45,
                    paddingAngle: 2,
                    cornerRadius: 3,
                    highlightScope: { fade: "global", highlight: "item" },
                  },
                ]}
                colors={[palette.palette.success.main, palette.palette.grey[400]]}
                height={220}
                slotProps={{ legend: { direction: "vertical" } }}
              />
            </ChartCard>
          </Box>

          <ChartCard
            title="Primary sequence length"
            description="Number of building blocks parsed per database entry"
          >
            <BarChart
              xAxis={[
                {
                  scaleType: "band",
                  data: stats.sequenceLengthBuckets.map((b) => b.label),
                  label: "Blocks in sequence",
                },
              ]}
              yAxis={[{ label: "Entries" }]}
              series={[
                {
                  data: stats.sequenceLengthBuckets.map((b) => b.count),
                  label: "Entries",
                  color: palette.palette.primary.main,
                },
              ]}
              height={280}
              hideLegend
              grid={{ horizontal: true }}
            />
          </ChartCard>

          <ChartCard
            title="Most common building blocks"
            description="Top identified monomers and tailoring events across the database"
          >
            <BarChart
              layout="horizontal"
              yAxis={[
                {
                  scaleType: "band",
                  data: [...stats.topBlocks].reverse().map((b) => b.label),
                },
              ]}
              xAxis={[{ label: "Occurrences" }]}
              series={[
                {
                  data: [...stats.topBlocks].reverse().map((b) => b.count),
                  label: "Occurrences",
                  color: palette.palette.primary.main,
                },
              ]}
              height={360}
              hideLegend
              grid={{ vertical: true }}
              margin={{ left: 110 }}
            />
          </ChartCard>
        </>
      )}
    </Box>
  );
};
