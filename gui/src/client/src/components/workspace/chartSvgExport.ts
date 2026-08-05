// Exports an on-screen custom SVG chart (violin, PMI shape triangle, ...) as a
// self-contained, print-quality SVG. Rather than re-deriving each chart's geometry,
// this clones the live rendered <svg> (so it can never drift from what's on screen)
// and patches only the handful of theme-dependent colors that would otherwise be
// unreadable outside the app -- gridlines/axis text sized for the current theme, and
// dot rings cut from the current surface color -- to fixed, white-background-friendly
// values. Mirrors alignmentSvgExport.ts's same print-quality convention for the
// alignment views. Any chart that wants this just tags its gridlines/axis-text/rings
// with the matching data-export-role attribute -- see ViolinChart/ShapeChart.

import { downloadSvg } from "./alignmentSvgExport";

const EXPORT_BACKGROUND = "#ffffff";
const EXPORT_GRID_COLOR = "#e0e0e0";
const EXPORT_AXIS_TEXT_COLOR = "#555555";
const EXPORT_RING_COLOR = "#ffffff";

const SVG_NS = "http://www.w3.org/2000/svg";

export function exportChartSvg(svgEl: SVGSVGElement, filename: string): void {
  const clone = svgEl.cloneNode(true) as SVGSVGElement;
  clone.setAttribute("xmlns", SVG_NS);

  const width = clone.getAttribute("width") ?? "0";
  const height = clone.getAttribute("height") ?? "0";

  const background = document.createElementNS(SVG_NS, "rect");
  background.setAttribute("x", "0");
  background.setAttribute("y", "0");
  background.setAttribute("width", width);
  background.setAttribute("height", height);
  background.setAttribute("fill", EXPORT_BACKGROUND);
  clone.insertBefore(background, clone.firstChild);

  clone.querySelectorAll('[data-export-role="gridline"]').forEach((el) => {
    el.setAttribute("stroke", EXPORT_GRID_COLOR);
  });
  clone.querySelectorAll('[data-export-role="axis-text"]').forEach((el) => {
    el.setAttribute("fill", EXPORT_AXIS_TEXT_COLOR);
  });
  clone.querySelectorAll('[data-export-role="dot-ring"]').forEach((el) => {
    el.setAttribute("stroke", EXPORT_RING_COLOR);
  });

  const markup = new XMLSerializer().serializeToString(clone);
  downloadSvg(markup, filename);
}
