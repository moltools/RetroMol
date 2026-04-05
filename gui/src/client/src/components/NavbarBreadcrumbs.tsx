import React from "react";
import { styled } from "@mui/material/styles";
import Breadcrumbs, { breadcrumbsClasses } from "@mui/material/Breadcrumbs";
import Typography from "@mui/material/Typography";
import NavigateNextRoundedIcon from "@mui/icons-material/NavigateNextRounded";
import { useLocation, Link } from "react-router-dom";

const StyledBreadcrumbs = styled(Breadcrumbs)(({ theme }) => ({
  margin: theme.spacing(1, 0),
  [`& .${breadcrumbsClasses.separator}`]: {
    color: theme.vars ? theme.vars.palette.action.disabled : theme.palette.action.disabled,
    margin: 1,
  },
  [`& .${breadcrumbsClasses.ol}`]: {
    alignItems: "center",
  },
}))

export const NavbarBreadcrumbs: React.FC = () => {
  const location = useLocation();
  const pathnames = location.pathname.split("/").filter(Boolean);

  return (
    <StyledBreadcrumbs aria-label="breadcrumb" separator={<NavigateNextRoundedIcon fontSize="small" />}>
      <Link to="/" style={{ textDecoration: "none", color: "inherit" }}>
        <Typography variant="body1">Home</Typography>
      </Link>
      {pathnames.map((value, index) => {
        const to = `/${pathnames.slice(0, index + 1).join("/")}`;
        const isLast = index === pathnames.length - 1;
        const label = value.charAt(0).toUpperCase() + value.slice(1);
        return isLast ? (
          <Typography key={to} variant="body1" sx={{ color: "text.primary", fontWeight: 600 }}>
            {label}
          </Typography>
        ) : (
          <Link key={to} to={to} style={{ textDecoration: "none", color: "inherit" }}>
            <Typography variant="body1">{label}</Typography>
          </Link>
        );
      })}
    </StyledBreadcrumbs>
  )
}