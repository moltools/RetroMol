import React from "react";
import Backdrop from "@mui/material/Backdrop";
import CircularProgress from "@mui/material/CircularProgress";

interface OverlayContextProps {
  showOverlay: () => void;
  hideOverlay: () => void;
}

// Create the context with an undefined initial value
const OverlayContext = React.createContext<OverlayContextProps | undefined>(undefined);

// Custom hook to use the overlay context
export const useOverlay = (): OverlayContextProps => {
  const context = React.useContext(OverlayContext);
  if (!context) {
    throw new Error("useOverlay must be used within an OverlayProvider");
  }
  return context;
}

interface OverlayProviderProps {
  children: React.ReactNode;
}

// The provider component wraps the app and renders the overlay
export const OverlayProvider: React.FC<OverlayProviderProps> = ({ children }) => {
  const [open, setOpen] = React.useState<boolean>(false);

  const showOverlay = React.useCallback(() => setOpen(true), []);
  const hideOverlay = React.useCallback(() => setOpen(false), []);

  return (
    <OverlayContext.Provider value={{ showOverlay, hideOverlay }}>
      {children}
      <Backdrop
        open={open}
        sx={{
          position: "fixed", // cover the entire viewport
          top: 0,
          left: 0,
          width: "100vw",
          height: "100vh",
          zIndex: (theme) => theme.zIndex.modal + 1, // ensure it sits on top
          color: "#fff",
          backgroundColor: "rgba(0, 0, 0, 0.5)", // fuzzy, foggy look
        }}
      >
        <CircularProgress color="inherit" />
      </Backdrop>
    </OverlayContext.Provider>
  )
}