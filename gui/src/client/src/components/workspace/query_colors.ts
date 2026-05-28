export const canonicalMotifKey = (s: string): string => {
  return s.trim().replace(/\s+/g, "").replace(/\^[RS]/g, "");
};

function parseColor(color: string, alpha: number): string {
  if (color.startsWith("#")) {
    let hex = color.replace(/^#/, "");
    if (hex.length === 3) hex = hex.split("").map((c) => c + c).join("");
    const r = parseInt(hex.slice(0, 2), 16);
    const g = parseInt(hex.slice(2, 4), 16);
    const b = parseInt(hex.slice(4, 6), 16);
    return `rgba(${r}, ${g}, ${b}, ${alpha})`;
  }

  const hsl = color.match(
    /hsl\(\s*([\d.]+)(?:deg)?\s*,\s*([\d.]+)%\s*,\s*([\d.]+)%\s*\)/
  );
  if (hsl) {
    const h = hsl[1];
    const s = hsl[2];
    const l = hsl[3];
    return `hsla(${h}, ${s}%, ${l}%, ${alpha})`;
  }

  throw new Error(`Unsupported color format: ${color}`);
}

export const defaultMotifColorMap = (): Record<string, string> => {
  const newColorMap: Record<string, string> = {};

  const baseColors: Record<"A" | "B" | "C" | "D", string> = {
    A: "#e74c3c",
    B: "#27ae60",
    C: "#2980b9",
    D: "#f39c12",
  };

  for (const key of Object.keys(baseColors) as Array<keyof typeof baseColors>) {
    const color = baseColors[key];
    newColorMap[key] = color;

    for (let i = 1; i <= 20; i++) {
      const alpha = 1 - i / 20;
      const alphaRounded = Math.round(alpha * 1000) / 1000;
      newColorMap[`${key}${i}`] = parseColor(color, alphaRounded);
    }
  }

  return newColorMap;
};

// Build once per module-load
const MOTIF_COLOR_MAP = defaultMotifColorMap();

export const getMotifColor = (name: string): string | null => {
  const key = canonicalMotifKey(name);
  return MOTIF_COLOR_MAP[key] || null;
};

// --- svg normalization helpers (moved from component) ---

const toHex = (v: number) => v.toString(16).padStart(2, "0");

const hslToRgb = (h: number, s: number, l: number) => {
  const c = (1 - Math.abs(2 * l - 1)) * s;
  const x = c * (1 - Math.abs(((h / 60) % 2) - 1));
  const m = l - c / 2;
  const pick = (hp: number) =>
    hp < 60
      ? [c, x, 0]
      : hp < 120
        ? [x, c, 0]
        : hp < 180
          ? [0, c, x]
          : hp < 240
            ? [0, x, c]
            : hp < 300
              ? [x, 0, c]
              : [c, 0, x];
  const [r1, g1, b1] = pick(h);
  return [
    Math.round((r1 + m) * 255),
    Math.round((g1 + m) * 255),
    Math.round((b1 + m) * 255),
  ];
};

export const normalizeColor = (raw: string | undefined | null) => {
  if (!raw) return "#f5f5f5";
  const c = raw.trim();

  if (/^#([0-9a-f]{3}|[0-9a-f]{6}|[0-9a-f]{8})$/i.test(c)) return c;

  const rgba = c.match(
    /^rgba?\(\s*([\d.]+)\s*,\s*([\d.]+)\s*,\s*([\d.]+)(?:\s*,\s*([\d.]+))?\s*\)$/i
  );
  if (rgba) {
    const [r, g, b, aRaw] = rgba.slice(1).map(Number);
    const a = isNaN(aRaw) ? 1 : Math.max(0, Math.min(1, aRaw));
    const blend = (v: number) => Math.round((1 - a) * 255 + a * v);
    return `#${toHex(blend(r))}${toHex(blend(g))}${toHex(blend(b))}`;
  }

  const hsla = c.match(
    /^hsla?\(\s*([\d.]+)(?:deg)?\s*,\s*([\d.]+)%\s*,\s*([\d.]+)%(?:\s*,\s*([\d.]+))?\s*\)$/i
  );
  if (hsla) {
    const h = Number(hsla[1]);
    const s = Number(hsla[2]) / 100;
    const l = Number(hsla[3]) / 100;
    const a =
      hsla[4] === undefined ? 1 : Math.max(0, Math.min(1, Number(hsla[4])));
    const [r, g, b] = hslToRgb(h, s, l);
    const blend = (v: number) => Math.round((1 - a) * 255 + a * v);
    return `#${toHex(blend(r))}${toHex(blend(g))}${toHex(blend(b))}`;
  }

  return "#f5f5f5";
};