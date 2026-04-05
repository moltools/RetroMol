import Box from "@mui/material/Box";
import Container from "@mui/material/Container";
import IconButton from "@mui/material/IconButton";
import Link from "@mui/material/Link";
import Stack from "@mui/material/Stack";
import Typography from "@mui/material/Typography";
import ArticleIcon from "@mui/icons-material/Article";
import GitHubIcon from "@mui/icons-material/GitHub";

function Copyright() {
  return (
    <Typography variant="body2" sx={{ color: "text.secondary", mt: 1 }}>
      {"Copyright © "}
      <Link
        color="text.secondary"
        onClick={
          () => {
              window.open(
                "https://github.com/moltools/RetroMol/issues",
                "_blank",
                "noopener,noreferrer"
              );
            }
        }
      >
        RetroMol
      </Link>
      &nbsp;
      {new Date().getFullYear()}
    </Typography>
  )
}

function FooterButtons() {
  return (
    <Stack
      direction="row"
      spacing={1}
      useFlexGap
      sx={{ justifyContent: "left", color: "text.secondary" }}
    >
      <IconButton
        color="inherit"
        size="small"
        aria-label="GitHub"
        sx={{ alignSelf: "center" }}
        onClick={
          () => {
              window.open(
                "https://github.com/moltools/RetroMol/issues",
                "_blank",
                "noopener,noreferrer"
              );
            }
        }
      >
        <GitHubIcon />
      </IconButton>
      <IconButton
        color="inherit"
        size="small"
        href=""
        aria-label="Article"
        sx={{ alignSelf: "center" }}
        disabled
      >
        <ArticleIcon />
      </IconButton>
    </Stack>
  )
}

export default function Footer() {
  return (
    <Container
      sx={{
        display: "flex",
        flexDirection: "column",
        alignItems: "center",
        py: { xs: 5 },
        textAlign: { sm: "center", md: "left" },
      }}
    >
      <Box sx={{ display: "flex", justifyContent: "space-between", width: "100%" }}>
        <Copyright />
        <FooterButtons />
      </Box>
    </Container>
  )
}