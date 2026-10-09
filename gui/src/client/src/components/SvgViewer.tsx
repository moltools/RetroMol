import React from "react";
import DOMPurify from "dompurify";
import Box from "@mui/material/Box";
import Paper from "@mui/material/Paper";
import Slider from "@mui/material/Slider";
import Stack from "@mui/material/Stack";
import Tooltip from "@mui/material/Tooltip";
import Typography from "@mui/material/Typography";
import AddIcon from "@mui/icons-material/Add";
import RemoveIcon from "@mui/icons-material/Remove";
import FitScreenIcon from "@mui/icons-material/FitScreen";
import DownloadIcon from "@mui/icons-material/Download";

export interface SvgViewerProps {
  svg: string; // raw SVG string
  initialZoom?: number; // initial zoom level (default is 1)
  minZoom?: number;
  maxZoom?: number;
  zoomStep?: number;
  onZoomChange?: (zoom: number) => void;
  onElementClick?: (info: {
    targetId: string | null;
    clientX: number;
    clientY: number;
  }) => void;
  height?: number | string; // optional fixed height
  downloadFileName?: string;
}

export const SvgViewer: React.FC<SvgViewerProps> = ({
  svg,
  initialZoom = 1,
  minZoom = 0.25,
  maxZoom = 4,
  zoomStep = 0.25,
  onZoomChange,
  onElementClick,
  height = 400,
  downloadFileName = "image.svg",
}) => {
  const [zoom, setZoom] = React.useState(initialZoom);
  const [isPanning, setIsPanning] = React.useState(false);
  const [pan, setPan] = React.useState({ x: 0, y: 0 });
  const [panStart, setPanStart] = React.useState<{ x: number; y: number } | null>(null);
  const [panOrigin, setPanOrigin] = React.useState<{ x: number; y: number }>({ x: 0, y: 0 });

  const containerRef = React.useRef<HTMLDivElement | null>(null);

  // Server-provided SVG markup is untrusted input (derived from user-uploaded
  // compounds/gene clusters) and gets injected raw via dangerouslySetInnerHTML,
  // so it must be sanitized first to prevent stored/reflected XSS.
  const sanitizedSvg = React.useMemo(
    () => DOMPurify.sanitize(svg, { USE_PROFILES: { svg: true, svgFilters: true } }),
    [svg]
  );

  const handleFit = React.useCallback(() => {
    // naive fit: reset zoom and pan
    setZoom(initialZoom);
    setPan({ x: 0, y: 0 });
    setPanOrigin({ x: 0, y: 0 });
  }, [initialZoom]);

  // Reset zoom and pan whenever SVG changes
  React.useEffect(() => {
    handleFit();
  }, [svg, handleFit]);

  // Notify parent about zoom changes
  React.useEffect(() => {
    onZoomChange?.(zoom);
  }, [zoom, onZoomChange])

  const clampZoom = (z: number) => Math.min(maxZoom, Math.max(minZoom, z));

  const handleZoomIn = () => {
    setZoom((z) => clampZoom(z + zoomStep));
  }

  const handleZoomOut = () => {
    setZoom((z) => clampZoom(z - zoomStep));
  }

  const handleZoomSliderChange = (event: Event, value: number | number[]) => {
    if (typeof value === "number") {
      setZoom(value);
    }
  }

  const handleMouseDown = (e: React.MouseEvent<HTMLDivElement>) => {
    if (e.button !== 0) return;
    setIsPanning(true);
    setPanStart({ x: e.clientX, y: e.clientY });
    setPanOrigin({ x: pan.x, y: pan.y });
  }

  const handleMouseMove = (e: React.MouseEvent<HTMLDivElement>) => {
    if (!isPanning || !panStart) return;
    const dx = e.clientX - panStart.x;
    const dy = e.clientY - panStart.y;
    setPan({ x: panOrigin.x + dx, y: panOrigin.y + dy });
  }

  const endPan = () => {
    setIsPanning(false);
    setPanStart(null);
  }

  const handleClick = (e: React.MouseEvent<HTMLDivElement>) => {
    if (!onElementClick) return;

    // Do not treat clicks that were actual pans as element clicks
    if (isPanning || panStart !== null) return;

    const target = e.target as HTMLElement | null;
    const targetId = target?.id ?? null;

    onElementClick({
      targetId,
      clientX: e.clientX,
      clientY: e.clientY,
    });
  }

  const handleDownload = () => {
    // Download the sanitized markup too — a raw SVG opened directly in a
    // browser tab can still execute embedded <script> content.
    const blob = new Blob([sanitizedSvg], { type: "image/svg+xml;charset=utf-8" });
    const url = URL.createObjectURL(blob);
    const link = document.createElement("a");
    link.href = url;
    link.download = downloadFileName;
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
    URL.revokeObjectURL(url);
  }

  return (
    <Paper
      elevation={3}
      sx={{
        display: "flex",
        flexDirection: "column",
        width: "100%",
        height,
        overflow: "hidden",
      }}
    >

      {/* Toolbar */}
      <Box
        sx={{
          px: 1.5,
          py: 1,
          borderBottom: (theme) => `1px solid ${theme.palette.divider}`,
        }}
      >
        <Stack direction="row" alignItems="center" spacing={1}>
          <Tooltip title="Zoom out" arrow>
            <RemoveIcon
              sx={{ cursor: "pointer" }}
              onClick={() => {
                // prevent zooming out beyond minZoom
                if (zoom <= minZoom) { return; }
                handleZoomOut();
              }}
            />
          </Tooltip>
          <Slider
            value={zoom}
            min={minZoom}
            max={maxZoom}
            step={zoomStep}
            onChange={handleZoomSliderChange}
            sx={{ width: 120 }}
          />
          <Tooltip title="Zoom in" arrow>
            <AddIcon
              sx={{ cursor: "pointer" }}
              onClick={() => {
                // prevent zooming in beyond maxZoom
                if (zoom >= maxZoom) { return; }
                handleZoomIn();
              }}
            />
          </Tooltip>
          <Tooltip title="Fit to screen" arrow>
            <FitScreenIcon
              sx={{ cursor: "pointer" }}
              onClick={handleFit}
            />
          </Tooltip>
          <Tooltip title="Download SVG" arrow>
            <DownloadIcon
              sx={{ cursor: "pointer" }}
              onClick={handleDownload}
            />
          </Tooltip>
          <Box flexGrow={1} />
          <Typography variant="caption" color="text.secondary">
            {Math.round(zoom * 100)}%
          </Typography>
        </Stack>
      </Box>

      {/* Viewport */}
      <Box
        ref={containerRef}
        onMouseDown={handleMouseDown}
        onMouseMove={handleMouseMove}
        onMouseUp={endPan}
        onMouseLeave={endPan}
        onClick={handleClick}
        sx={{
          flexGrow: 1,
          position: "relative",
          overflow: "hidden",
          cursor: isPanning ? "grabbing" : "grab",
          backgroundColor: (theme) => theme.palette.background.default,
          userSelect: isPanning ? "none" : "auto",
        }}
      >
        <Box
          sx={{
            transform: `translate(${pan.x}px, ${pan.y}px) scale(${zoom})`,
            transformOrigin: "0 0",
            width: "fit-content",
            height: "fit-content",
            // Required so clicks hit SVG and not the wrapper
            "& svg": {
              display: "block",
            }
          }}
          dangerouslySetInnerHTML={{ __html: sanitizedSvg }}
        />
      </Box>

    </Paper>
  )
}