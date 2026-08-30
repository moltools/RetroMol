import Box from "@mui/material/Box";
import Container from "@mui/material/Container";
import Stack from "@mui/material/Stack";
import Typography from "@mui/material/Typography";
import { ServerUptime } from "./ServerUptime";

export default function Hero() {
  return (
    <Box id="hero" sx={{ width: "100%" }}>
      <Container
        id="retrieve-session"
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
          <Box
            component="img"
            src={`${process.env.PUBLIC_URL}/logo.svg`}
            alt="RetroMol logo"
            sx={{
              width: 300,
              height: "auto",
              position: "relative",
              zIndex: 1,

              // Stacking drop shadows to mimic an outline that follows the SVG shape
              filter: `
                drop-shadow(2px 2px 0 white)
                drop-shadow(-2px 2px 0 white)
                drop-shadow(2px -2px 0 white)
                drop-shadow(-2px -2px 0 white)
              `,
            }}
          />
          <Typography
            variant="h1"
            sx={{
              display: "flex",
              flexDirection: { xs: "column", sm: "row" },
              alignItems: "center",
              fontSize: "clamp(3rem, 10vw, 3.5rem)",
            }}
          >
            Discover&nbsp;
            <Typography
              component="span"
              variant="h1"
              sx={(theme) => ({
                fontSize: "inherit",
                color: "primary.main",
                ...theme.applyStyles("dark", {
                  color: "primary.light",
                }),
              })}
            >
              RetroMol
            </Typography>
          </Typography>
          <Typography
            sx={{
              textAlign: "center",
              color: "text.secondary",
              width: { sm: "100%", md: "80%" },
            }}
          >
            Perform cross-modal retrieval between natural product compounds and
            BGCs.
          </Typography>
          <ServerUptime />
        </Stack>
      </Container>
    </Box>
  )
}