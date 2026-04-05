export function createCookie(name: string, value: string, days?: number, path: string = "/"): void {
  let cookieStr = `${encodeURIComponent(name)}=${encodeURIComponent(value)}`;
  if (typeof days === "number") {
    const expires = new Date();
    // Set the expiration date for the cookie; here we set it to "days" days from now
    expires.setTime(expires.getTime() + days * 24 * 60 * 60 * 1000);
    cookieStr += `; expires=${expires.toUTCString()}`;
  }
  cookieStr += `; path=${path}`;
  document.cookie = cookieStr;
}

export function getCookie(name: string): string | null {
  const match = document.cookie
    .split("; ")
    .map(pair => pair.split("="))
    .find(([key]) => key === name);
  return match ? match[1] : null;
}

export function deleteCookie(name: string, path: string = "/"): void {
  document.cookie =
    `${encodeURIComponent(name)}=; ` +
    `expires=Thu, 01 Jan 1970 00:00:00 GMT; ` +
    `path=${path};`;
}

function parseColor(color: string, alpha: number): string {
  // HEX case: "#RGB" or "#RRGGBB"
  if (color.startsWith("#")) {
    let hex = color.replace(/^#/, "");
    // expand shorthand (#abc → aabbcc)
    if (hex.length === 3) {
      hex = hex.split("").map(c => c + c).join("");
    }
    // parse r, g, b
    const r = parseInt(hex.slice(0, 2), 16);
    const g = parseInt(hex.slice(2, 4), 16);
    const b = parseInt(hex.slice(4, 6), 16);
    return `rgba(${r}, ${g}, ${b}, ${alpha})`;
  };

  // HSL case: "hsl(h, s%, l%)"
  const hsl = color.match(
    /hsl\(\s*([\d.]+)(?:deg)?\s*,\s*([\d.]+)%\s*,\s*([\d.]+)%\s*\)/
  );
  if (hsl) {
    const h = hsl[1];
    const s = hsl[2];
    const l = hsl[3];
    return `hsla(${h}, ${s}%, ${l}%, ${alpha})`;
  };

  throw new Error(`Unsupported color format: ${color}`);
};

export const defaultMotifColorMap = (): Record<string, string> => {
  const newColorMap: Record<string, string> = {};

  const baseColors: Record<"A"|"B"|"C"|"D", string> = {
    A: "#e74c3c ", // red
    B: "#27ae60", // green
    C: "#2980b9", // blue
    D: "#f39c12 ", // orange
  };

  for (const key of Object.keys(baseColors) as Array<keyof typeof baseColors>) {
    const color = baseColors[key];
    // plain (opaque) base
    newColorMap[key] = color;

    // numbered variants 1→15 → alpha = 1/15…15/15
    for (let i = 1; i <= 15; i++) {
      const alpha = 1 - (i / 15);
      const alphaRounded = Math.round(alpha * 1000) / 1000;
      newColorMap[`${key}${i}`] = parseColor(color, alphaRounded);
    }
  }

  return newColorMap;
};