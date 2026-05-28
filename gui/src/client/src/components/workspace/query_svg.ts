import { MsaRow, Sequence } from "./query_types";
import { getMotifColor, normalizeColor } from "./query_colors";
import { isPolyketideMotif } from "./query_motif";

export const escapeSvgText = (value: string) =>
  value
    .replace(/&/g, "&amp;")
    .replace(/</g, "&lt;")
    .replace(/>/g, "&gt;")
    .replace(/"/g, "&quot;")
    .replace(/'/g, "&#39;");

type RichLabel =
  | { isRich: false; text: string }
  | { isRich: true; parts: string[] };

const svgMotifLabel = (
  rawName: string | null,
  toDisplayName: (name: string | null) => string | null
): RichLabel => {
  const raw = (rawName || "").trim();
  if (!raw) return { isRich: false, text: "UNK" };

  if (isPolyketideMotif(raw)) {
    const parts = raw.split(/(\^[RS])/g).filter(Boolean);
    return { isRich: true, parts };
  }

  return { isRich: false, text: toDisplayName(raw) || "X" };
};

const renderSvgTspansForPolyketide = (parts: string[]) => {
  return parts
    .map((p) => {
      if (p === "^R" || p === "^S") {
        const v = escapeSvgText(p.slice(1));
        return `<tspan baseline-shift="super" font-size="5.2">${v}</tspan>`;
      }
      return `<tspan>${escapeSvgText(p)}</tspan>`;
    })
    .join("");
};

export const sequenceLength = (seqs: Sequence[]) =>
  seqs.reduce((sum, seq) => sum + seq.sequence.length, 0);

type BuildSvgArgs = {
  msa: MsaRow[];
  msaLength: number;
  toDisplayName: (name: string | null) => string | null;
};

export const buildMsaSvg = ({ msa, msaLength, toDisplayName }: BuildSvgArgs) => {
  const motifPx = 19.1955;
  const rowChipH = 13.8131;

  const rowGap = 5;
  const subHeaderH = 8;
  const rowBlockH = subHeaderH + rowChipH;

  const padding = 12;

  const labelW = 120;
  const scoreW = 42;
  const leftW = labelW + scoreW;

  const svgWidth = padding * 2 + leftW + motifPx * msaLength;
  const svgHeight =
    padding * 2 + msa.length * rowBlockH + Math.max(0, msa.length - 1) * rowGap;

  const lineSpan = Math.max(0, msaLength - 1) * motifPx;

  const fmt = (v: number | null, digits: number) =>
    v == null || Number.isNaN(v) ? "" : v.toFixed(digits);

  const rowsSvg = msa
    .map((row, rIdx) => {
      const yTop = padding + rIdx * (rowBlockH + rowGap);
      const yHeader = yTop;
      const yChips = yTop + subHeaderH;

      const rowName = escapeSvgText(row.name || row.id || "row");

      const alignText = escapeSvgText(fmt(row.alignment_score, 2));
      const cosineText = escapeSvgText(fmt(row.cosine_score, 2));
      const matchText = escapeSvgText(fmt(row.match_score, 2));

      const xLeft = padding;
      const xName = xLeft + 6;
      const xScoresRight = xLeft + leftW - 6;
      const xMotifs = xLeft + leftW;

      const lineY = yChips + rowChipH / 2;
      const lineX1 = xMotifs;
      const lineX2 = xMotifs + lineSpan + motifPx;

      let col = 0;

      const subseqHeaders = row.sequence
        .map((subseq) => {
          const startCol = col;
          const len = subseq.sequence.length;
          col += len;

          const allGaps = subseq.sequence.every((it) => it.isGap);
          if (allGaps) return "";

          const x = xMotifs + startCol * motifPx;
          const w = len * motifPx;
          const title = escapeSvgText(subseq.name || subseq.id || "");

          const h = subHeaderH + rowChipH / 2;
          const r = 4;

          const headerPath = `
            M ${x} ${yHeader + h}
            L ${x} ${yHeader + r}
            Q ${x} ${yHeader} ${x + r} ${yHeader}
            L ${x + w - r} ${yHeader}
            Q ${x + w} ${yHeader} ${x + w} ${yHeader + r}
            L ${x + w} ${yHeader + h}
            L ${x} ${yHeader + h}
          `;

          return `
            <g>
              <path
                d="${headerPath}"
                fill="#e0e0e0"
                stroke="#bdbdbd"
                stroke-width="0.8"
              />
              <text x="${x + 3}" y="${yHeader + subHeaderH / 2 + 1}"
                font-family="Helvetica" font-size="4" fill="#111"
                dominant-baseline="middle">${title}</text>
            </g>
          `;
        })
        .join("");

      col = 0;

      const cells = row.sequence
        .flatMap((subseq) =>
          subseq.sequence.map((motif) => {
            const x = xMotifs + col * motifPx;
            col += 1;

            if (motif.isGap) return "";

            const fill = normalizeColor(getMotifColor(motif.name || "") || "#ffffff");
            const label = svgMotifLabel(motif.name, toDisplayName);
            const paddingText = 3;

            return `
              <g>
                <rect x="${x}" y="${yChips}" width="${motifPx}" height="${rowChipH}"
                  rx="4" ry="4" fill="${fill}" stroke="#000000" stroke-width="0.9" />
                <text x="${x + motifPx / 2}" y="${yChips + rowChipH / 2 + paddingText}"
                  font-family="Helvetica" font-size="7.2" font-weight="600" fill="#000"
                  dominant-baseline="middle" text-anchor="middle">
                  ${
                    label.isRich
                      ? renderSvgTspansForPolyketide(label.parts)
                      : escapeSvgText(label.text)
                  }
                </text>
              </g>
            `;
          })
        )
        .join("");

      return `
        <g>
          <text x="${xName}" y="${yChips + rowChipH / 2 + 3}"
            font-family="Helvetica" font-size="8.2" font-weight="600"
            dominant-baseline="middle">${rowName}</text>

          <text x="${xScoresRight}" y="${yChips + 3.5}"
            font-family="monospace" font-size="7.0" fill="#111"
            text-anchor="end">${alignText}</text>
          <text x="${xScoresRight}" y="${yChips + 9.0}"
            font-family="monospace" font-size="7.0" fill="#111"
            text-anchor="end">${cosineText}</text>
          <text x="${xScoresRight}" y="${yChips + 14.5}"
            font-family="monospace" font-size="7.0" fill="#111"
            text-anchor="end">${matchText}</text>

          ${
            lineSpan > 0
              ? `<line x1="${lineX1}" x2="${lineX2}" y1="${lineY}" y2="${lineY}"
                  stroke="#9e9e9e" stroke-width="0.9" />`
              : ""
          }

          ${subseqHeaders}
          ${cells}
        </g>
      `;
    })
    .join("");

  return `
    <svg xmlns="http://www.w3.org/2000/svg" width="${svgWidth}" height="${svgHeight}">
      <rect x="0" y="0" width="${svgWidth}" height="${svgHeight}" fill="#ffffff"/>
      ${rowsSvg}
    </svg>`;
};