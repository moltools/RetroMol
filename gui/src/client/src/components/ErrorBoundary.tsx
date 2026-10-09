import React from "react";
import Alert from "@mui/material/Alert";
import AlertTitle from "@mui/material/AlertTitle";
import Button from "@mui/material/Button";
import Box from "@mui/material/Box";

interface ErrorBoundaryProps {
  children: React.ReactNode;
  /** Rendered instead of the default alert when a child throws. */
  fallback?: React.ReactNode;
  /** Short label for the default fallback, e.g. "molecule structure". */
  what?: string;
}

interface ErrorBoundaryState {
  error: Error | null;
}

// React error boundaries must be class components — there is no hook equivalent.
// Catches render/lifecycle errors in the subtree (e.g. a third-party drawing
// library choking on a malformed user-uploaded SMILES/GenBank file) so a single
// bad input can't white-screen the whole app.
export class ErrorBoundary extends React.Component<ErrorBoundaryProps, ErrorBoundaryState> {
  state: ErrorBoundaryState = { error: null };

  static getDerivedStateFromError(error: Error): ErrorBoundaryState {
    return { error };
  }

  componentDidCatch(error: Error, info: React.ErrorInfo) {
    console.error("ErrorBoundary caught an error:", error, info.componentStack);
  }

  private reset = () => this.setState({ error: null });

  render() {
    if (this.state.error) {
      if (this.props.fallback !== undefined) return this.props.fallback;

      return (
        <Box sx={{ p: 2 }}>
          <Alert severity="error" action={
            <Button color="inherit" size="small" onClick={this.reset}>
              Retry
            </Button>
          }>
            <AlertTitle>Something went wrong{this.props.what ? ` rendering the ${this.props.what}` : ""}</AlertTitle>
            {this.state.error.message || "An unexpected error occurred."}
          </Alert>
        </Box>
      );
    }

    return this.props.children;
  }
}
