import React from "react";
import ReactDOM from "react-dom/client";
import { StyledEngineProvider } from "@mui/material/styles";
import type {} from "@mui/material/themeCssVarsAugmentation";
import { BrowserRouter, Route, Routes, useLocation } from "react-router-dom";
import { AnimatePresence, motion } from "framer-motion";

// Import pages
import Home from "./pages/Home";
import Dashboard from "./pages/Dashboard";
import NotFound from "./pages/NotFound";

// Generate and store a random deviceId if it doesn"t exist already.
const deviceIdKey = "deviceId";
if (!localStorage.getItem(deviceIdKey)) {
  const newDeviceId = Math.random().toString(36).substring(2, 15); // random string
  localStorage.setItem(deviceIdKey, newDeviceId);
};

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
      <Routes location={location} key={getPageKey(location.pathname)}>
        <Route path="/" element={<Page><Home /></Page>} />
        <Route path="/dashboard/*" element={<Page><Dashboard /></Page>} />
        <Route path="/maxsessions" element={<Page><NotFound title="Maximum sessions reached" subtitle="The maximum number of active sessions has been reached. Please try again later." /></Page>} />
        <Route path="/oops" element={<Page><NotFound title="Something went wrong" subtitle="Something catastrophic happened. Please reach out to us on GitHub." /></Page>} />
        <Route path="/notfound" element={<Page><NotFound title="Session not found" subtitle="The session you are looking for does not exist." /></Page>} />
        <Route path="/*" element={<Page><NotFound /></Page>} />
      </Routes>
    </AnimatePresence>
  );
};

function App() {
  return (
    <StyledEngineProvider injectFirst>
      <BrowserRouter>
        <AnimatedRoutes />
      </BrowserRouter>
    </StyledEngineProvider>
  );
};

ReactDOM.createRoot(document.querySelector("#root") as HTMLElement).render(
  <React.StrictMode>
    <App />
  </React.StrictMode>
);