import React from "react";
import Box from "@mui/material/Box";
import Button from "@mui/material/Button";
import Container from "@mui/material/Container";
import CssBaseline from "@mui/material/CssBaseline";
import Stack from "@mui/material/Stack";
import Typography from "@mui/material/Typography";
import AppTheme from "../theme/AppTheme";
import Footer from "../components/Footer";

// Helper function to pad numbers with a leading zero if needed
const pad = (num: number): string => num.toString().padStart(2, "0")

// Helper function to format seconds into days:hours:minutes:seconds
const formatUptime = (uptime: number): string => {
  const days = Math.floor(uptime / (24 * 3600));
  const hours = Math.floor((uptime % (24 * 3600)) / 3600);
  const minutes = Math.floor((uptime % 3600) / 60);
  const seconds = uptime % 60;
  return `${days} day${days !== 1 ? "s" : ""} and ${pad(hours)}:${pad(minutes)}:${pad(seconds)} hours`;
}

const Hero = ({ title, subtitle }: { title?: string; subtitle?: string }) => {
  const [uptime, setUptime] = React.useState(0);

  // Fetch uptime from api/startup endpoint
  React.useEffect(() => {
    const fetchUptime = async () => {
      try {
        const response = await fetch("/api/startup");
        if (!response.ok) {
          throw new Error("network response was not ok");
        }
        const data = await response.json();
        setUptime(data.uptime);
      } catch (err) {
        console.error("error fetching uptime:", err);
      }
    };
    fetchUptime();
  }, [])

  return (
    <Box id="hero" sx={{width: "100%"}}>
      <Container
        sx={{
          display: "flex",
          flexDirection: "column",
          alignItems: "center",
          pt: { xs: 14, sm: 20 },
          pb: { xs: 8, sm: 12 },
        }}
      >
        <Stack
          spacing={2}
          useFlexGap
          sx={{ alignItems: "center", width: { xs: "100%", sm: "70%" } }}
        >
          <Typography
            variant="h1"
            sx={{
              display: "flex",
              flexDirection: { xs: "column", sm: "row" },
              alignItems: "center",
              fontSize: "clamp(3rem, 10vw, 3.5rem)",
            }}
          >
            {title || "Page not found"}
          </Typography>
          <Typography
            sx={{
              textAlign: "center",
              color: "text.secondary",
              width: { sm: "100%", md: "80%" },
            }}
          >
            {subtitle || "The page you are looking for does not exist."}
          </Typography>
          <Typography
            sx={{
              textAlign: "center",
              color: "text.secondary",
              fontWeight: 800,
              width: { sm: "100%", md: "80%" },
            }}
          >
            {`Server uptime: ${formatUptime(uptime)}`}
          </Typography>
          <Button
            variant="contained"
            size="large"
            sx={{ width: { xs: "100%", sm: "auto" } }}
            onClick={() => {
              window.location.href = "/";
            }}
          >
            Go to home
          </Button>
        </Stack>
      </Container>
    </Box>
  )
}

export default function NotFound(props: { title?: string; subtitle?: string; disableCustomTheme?: boolean }) {
  return (
    <AppTheme {...props}>
      <CssBaseline enableColorScheme />
      <Box sx={{ display: "flex", flexDirection: "column", minHeight: "100vh" }}>
        <Box sx={{ flexGrow: 1 }}>
          <Hero title={props.title} subtitle={props.subtitle} />
        </Box>
        <Footer />
      </Box>
    </AppTheme>
  )
}