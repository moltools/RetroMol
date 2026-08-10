// Builds a self-contained, print-quality SVG of an alignment (pairwise or MSA) for
// inclusion in a manuscript, and triggers a browser download of it. Deliberately
// constructs real <rect>/<text> elements rather than rasterizing the on-screen DOM
// (e.g. via html2canvas), so the output stays crisp at any zoom/print size and opens
// cleanly in vector editors like Illustrator/Inkscape.

const FONT_FAMILY = "Arial, Helvetica, sans-serif";
const NAME_FONT_SIZE = 12;
const LABEL_FONT_SIZE = 12;
const SCORE_FONT_SIZE = 11;
const ROW_HEIGHT = 24;
const ROW_GAP = 4;
const CELL_PADDING_X = 6;
const CELL_GAP_X = 3;
const MARGIN = 12;
const MIN_CELL_WIDTH = 22;

export type AlignmentSvgRow = {
  label: string;
  score?: string; // pre-formatted for display, e.g. "82.3%" — omit for rows with no score (e.g. the query)
  // Per-column Tanimoto similarity (0-1) against the query's unit in that column, mirroring
  // the on-screen shading; index-aligned with `sequence`. Omit entirely for a neutral fill
  // throughout (e.g. the query row).
  matchStrengths?: (number | null)[];
  sequence: (string | null)[]; // null = gap
};

// Mirrors the on-screen shading scale (see similarityColor in WorkspaceDiscovery.tsx),
// but computed as a flat hex since this file has no access to the live MUI theme.
const MATCH_SHADE_MIN_ALPHA = 0.12;
const MATCH_SHADE_MAX_ALPHA = 0.65;
const MATCH_COLOR_RGB = [25, 118, 210]; // MUI default primary.main (#1976d2)

function matchFillColor(similarity: number | null | undefined): string {
  if (similarity === null || similarity === undefined) return "#ffffff";
  const clamped = Math.max(0, Math.min(1, similarity));
  const a = MATCH_SHADE_MIN_ALPHA + clamped * (MATCH_SHADE_MAX_ALPHA - MATCH_SHADE_MIN_ALPHA);
  const [r, g, b] = MATCH_COLOR_RGB.map((channel) => Math.round(channel * a + 255 * (1 - a)));
  return `rgb(${r}, ${g}, ${b})`;
}

function escapeXml(value: string): string {
  return value
    .replace(/&/g, "&amp;")
    .replace(/</g, "&lt;")
    .replace(/>/g, "&gt;")
    .replace(/"/g, "&quot;")
    .replace(/'/g, "&apos;");
}

// Splits a name like "A2^R" or "C^E2" into plain/superscript runs, matching MotifName's convention.
function splitStereoMarkers(name: string): { text: string; sup: boolean }[] {
  return name
    .split(/(\^[SREZ])/g)
    .filter((part) => part !== "")
    .map((part) =>
      part === "^S" || part === "^R" || part === "^E" || part === "^Z"
        ? { text: part.slice(1), sup: true }
        : { text: part, sup: false }
    );
}

function measureName(ctx: CanvasRenderingContext2D, name: string, fontSize: number): number {
  let width = 0;
  for (const part of splitStereoMarkers(name)) {
    ctx.font = `${part.sup ? fontSize * 0.7 : fontSize}px ${FONT_FAMILY}`;
    width += ctx.measureText(part.text).width;
  }
  return width;
}

function nameToTspans(name: string, fontSize: number): string {
  return splitStereoMarkers(name)
    .map((part) =>
      part.sup
        ? `<tspan baseline-shift="super" font-size="${(fontSize * 0.7).toFixed(1)}">${escapeXml(part.text)}</tspan>`
        : `<tspan>${escapeXml(part.text)}</tspan>`
    )
    .join("");
}

export function buildAlignmentSvg(rows: AlignmentSvgRow[]): string {
  if (rows.length === 0) {
    return `<svg xmlns="http://www.w3.org/2000/svg" width="1" height="1"></svg>`;
  }

  const canvas = document.createElement("canvas");
  const ctx = canvas.getContext("2d");
  if (!ctx) {
    throw new Error("Could not get a 2D canvas context for SVG text layout");
  }

  const columnCount = Math.max(...rows.map((r) => r.sequence.length));
  const hasScores = rows.some((r) => r.score !== undefined);

  ctx.font = `bold ${LABEL_FONT_SIZE}px ${FONT_FAMILY}`;
  let labelColWidth = 0;
  for (const row of rows) {
    labelColWidth = Math.max(labelColWidth, ctx.measureText(row.label).width);
  }
  labelColWidth += CELL_PADDING_X * 2;

  let scoreColWidth = 0;
  if (hasScores) {
    ctx.font = `${SCORE_FONT_SIZE}px ${FONT_FAMILY}`;
    for (const row of rows) {
      if (row.score !== undefined) {
        scoreColWidth = Math.max(scoreColWidth, ctx.measureText(row.score).width);
      }
    }
    scoreColWidth += CELL_PADDING_X * 2;
  }

  const columnWidths: number[] = [];
  for (let col = 0; col < columnCount; col++) {
    let width = MIN_CELL_WIDTH;
    for (const row of rows) {
      const name = row.sequence[col] ?? null;
      if (name !== null) {
        width = Math.max(width, measureName(ctx, name, NAME_FONT_SIZE) + CELL_PADDING_X * 2);
      }
    }
    columnWidths.push(width);
  }

  const sequenceAreaWidth =
    columnWidths.reduce((sum, w) => sum + w, 0) + Math.max(columnWidths.length - 1, 0) * CELL_GAP_X;
  const totalWidth =
    MARGIN * 2 +
    labelColWidth +
    CELL_GAP_X +
    sequenceAreaWidth +
    (hasScores ? CELL_GAP_X + scoreColWidth : 0);
  const totalHeight = MARGIN * 2 + rows.length * ROW_HEIGHT + Math.max(rows.length - 1, 0) * ROW_GAP;

  const svg: string[] = [];
  svg.push(
    `<svg xmlns="http://www.w3.org/2000/svg" width="${Math.ceil(totalWidth)}" height="${Math.ceil(
      totalHeight
    )}" viewBox="0 0 ${Math.ceil(totalWidth)} ${Math.ceil(totalHeight)}" font-family="${FONT_FAMILY}">`
  );
  svg.push(`<rect x="0" y="0" width="${Math.ceil(totalWidth)}" height="${Math.ceil(totalHeight)}" fill="#ffffff" />`);

  rows.forEach((row, rowIdx) => {
    const y = MARGIN + rowIdx * (ROW_HEIGHT + ROW_GAP);
    const textY = y + ROW_HEIGHT / 2 + NAME_FONT_SIZE * 0.35;

    let x = MARGIN;
    svg.push(
      `<text x="${x}" y="${textY.toFixed(
        1
      )}" font-size="${LABEL_FONT_SIZE}" font-weight="700" fill="#111111">${escapeXml(row.label)}</text>`
    );
    x += labelColWidth + CELL_GAP_X;

    for (let col = 0; col < columnCount; col++) {
      const width = columnWidths[col];
      const name = row.sequence[col] ?? null;
      const isGap = name === null;

      svg.push(
        `<rect x="${x.toFixed(1)}" y="${y}" width="${width.toFixed(1)}" height="${ROW_HEIGHT}" rx="3" ` +
          `fill="${isGap ? "#f2f2f2" : matchFillColor(row.matchStrengths?.[col])}" stroke="${isGap ? "#dddddd" : "#333333"}" stroke-width="1" ` +
          `${isGap ? 'stroke-dasharray="2,2"' : ""} />`
      );

      if (!isGap) {
        const cx = x + width / 2;
        svg.push(
          `<text x="${cx.toFixed(1)}" y="${textY.toFixed(
            1
          )}" font-size="${NAME_FONT_SIZE}" text-anchor="middle" fill="#111111">${nameToTspans(
            name,
            NAME_FONT_SIZE
          )}</text>`
        );
      }

      x += width + CELL_GAP_X;
    }

    if (hasScores && row.score !== undefined) {
      svg.push(
        `<text x="${x.toFixed(1)}" y="${textY.toFixed(
          1
        )}" font-size="${SCORE_FONT_SIZE}" fill="#444444">${escapeXml(row.score)}</text>`
      );
    }
  });

  svg.push(`</svg>`);
  return svg.join("\n");
}

export function downloadSvg(svgMarkup: string, filename: string): void {
  const blob = new Blob([svgMarkup], { type: "image/svg+xml;charset=utf-8" });
  const url = URL.createObjectURL(blob);
  const link = document.createElement("a");
  link.href = url;
  link.download = filename;
  document.body.appendChild(link);
  link.click();
  document.body.removeChild(link);
  URL.revokeObjectURL(url);
}
