import React, { Suspense } from "react";
import ReactDOM from "react-dom/client";
import { StyledEngineProvider } from "@mui/material/styles";
import type {} from "@mui/material/themeCssVarsAugmentation";
import { BrowserRouter, Route, Routes, useLocation } from "react-router-dom";
import { AnimatePresence, motion } from "framer-motion";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import Box from "@mui/material/Box";
import CircularProgress from "@mui/material/CircularProgress";

import { ErrorBoundary } from "./components/ErrorBoundary";
import NotFound from "./pages/NotFound";

// Home and Dashboard are code-split: Dashboard alone pulls in the charting,
// molecule-drawing, and workspace UI, none of which the landing page needs.
const Home = React.lazy(() => import("./pages/Home"));
const Dashboard = React.lazy(() => import("./pages/Dashboard"));

const queryClient = new QueryClient({
  defaultOptions: {
    queries: {
      retry: 1,
      refetchOnWindowFocus: false,
    },
  },
});

const RouteFallback = () => (
  <Box sx={{ display: "flex", alignItems: "center", justifyContent: "center", height: "100vh" }}>
    <CircularProgress />
  </Box>
);

// Wraps pages with an animation effect
const Page = ({ children }: { children: React.ReactNode }) => {
  return (
    <motion.div
      initial={{ opacity: 0, x: -50 }}
      animate={{ opacity: 1, x: 0 }}
      exit={{ opacity: 0, x: 50 }}
      transition={{ duration: 0.5 }}
      style={{ height: "100%" }}
    >
      {children}
    </motion.div>
  );
};

// Derive a constant key for AnimatePresence so that navigating within /dashboard doesn’t re-trigger the animation
const getPageKey = (pathname: string): string => {
  if (pathname.startsWith("/dashboard")) {
    return "/dashboard";
  }
  return pathname;
};

const AnimatedRoutes = () => {
  const location = useLocation();

  return (
    <AnimatePresence mode="wait">
      <Suspense fallback={<RouteFallback />}>
        <Routes location={location} key={getPageKey(location.pathname)}>
          <Route path="/" element={<Page><Home /></Page>} />
          <Route path="/dashboard/*" element={<Page><Dashboard /></Page>} />
          <Route path="/maxsessions" element={<Page><NotFound title="Maximum sessions reached" subtitle="The maximum number of active sessions has been reached. Please try again later." /></Page>} />
          <Route path="/oops" element={<Page><NotFound title="Something went wrong" subtitle="Something catastrophic happened. Please reach out to us on GitHub." /></Page>} />
          <Route path="/notfound" element={<Page><NotFound title="Session not found" subtitle="The session you are looking for does not exist." /></Page>} />
          <Route path="/*" element={<Page><NotFound /></Page>} />
        </Routes>
      </Suspense>
    </AnimatePresence>
  );
};

function App() {
  return (
    <QueryClientProvider client={queryClient}>
      <StyledEngineProvider injectFirst>
        <BrowserRouter>
          <ErrorBoundary
            fallback={<NotFound title="Something went wrong" subtitle="Something catastrophic happened. Please reach out to us on GitHub." />}
          >
            <AnimatedRoutes />
          </ErrorBoundary>
        </BrowserRouter>
      </StyledEngineProvider>
    </QueryClientProvider>
  );
};

ReactDOM.createRoot(document.querySelector("#root") as HTMLElement).render(
  <React.StrictMode>
    <App />
  </React.StrictMode>
);
