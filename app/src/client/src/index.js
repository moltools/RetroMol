import React from "react";
import ReactDOM from "react-dom/client";
import { BrowserRouter, Routes, Route } from "react-router-dom";

import CineMol from "./pages/CineMol";
import RetroMol from "./pages/RetroMol";
import LandingPage from "./pages/LandingPage";
import NotFoundPage from "./pages/NotFoundPage";
import Toast from "./components/common/Toast";

import "./index.css";

function AppRoutes () {
    return (
        <div>
            <Routes>
                <Route path="/" element={<LandingPage />} />
                <Route path="/cinemol" element={<CineMol />} />
                <Route path="/retromol" element={<RetroMol />} />
                <Route path="/*" element={<NotFoundPage />} />
            </Routes>
        </div>
    );
};

function App () {
    return (
        <div>
            <BrowserRouter>
                <AppRoutes />
                <Toast />
            </BrowserRouter>
        </div>
    );
};

const root = ReactDOM.createRoot(document.getElementById("root"));

root.render(<App />);