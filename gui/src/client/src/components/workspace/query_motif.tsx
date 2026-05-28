import React from "react";

export const PROTECTED_NAME_TO_CODE: Record<string, string> = {
  ALANINE: "ALA",
  CYSTEINE: "CYS",
  ASPARTICACID: "ASP",
  GLUTAMICACID: "GLU",
  PHENYLALANINE: "PHE",
  GLYCINE: "GLY",
  HISTIDINE: "HIS",
  ISOLEUCINE: "ILE",
  LYSINE: "LYS",
  LEUCINE: "LEU",
  METHIONINE: "MET",
  ASPARAGINE: "ASN",
  PROLINE: "PRO",
  GLUTAMINE: "GLN",
  ARGININE: "ARG",
  SERINE: "SER",
  THREONINE: "THR",
  VALINE: "VAL",
  TRYPTOPHAN: "TRP",
  TYROSINE: "TYR",
};

export const renderChiralSuperscripts = (label: string) => {
  const parts = label.split(/(\^[RS])/g).filter(Boolean);

  return (
    <>
      {parts.map((p, i) => {
        if (p === "^R" || p === "^S") {
          return (
            <sup key={i} style={{ fontSize: "0.7em", lineHeight: 0 }}>
              {p.slice(1)}
            </sup>
          );
        }
        return <React.Fragment key={i}>{p}</React.Fragment>;
      })}
    </>
  );
};

export const isPolyketideMotif = (s: string | null | undefined) => {
  if (!s) return false;
  return /^[A-D](\^[RS])*(\d+)?(\^[RS])*$/i.test(s.trim());
};

export const makeToDisplayName = (protectedNameToCode: Record<string, string>) => {
  const norm = (s: string) => s.replace(/[^a-z0-9]/gi, "").toUpperCase();

  const prot = new Map<string, string>(
    Object.entries(protectedNameToCode).map(([k, v]) => [norm(k), norm(v)])
  );
  const reserved = new Set<string>(Array.from(prot.values()));
  const used = new Set<string>(reserved);
  const cache = new Map<string, string>();

  const candidates = (s: string) => {
    const out: string[] = [];
    if (s.length >= 3) {
      out.push(s.slice(0, 3));
      for (let i = 3; i < s.length; i++) out.push(s[0] + s[1] + s[i]);
      for (let i = 2; i < s.length; i++) out.push(s[0] + s[i - 1] + s[i]);
    }
    if (s.length >= 2) out.push(s.slice(0, 2));
    if (s.length >= 1) out.push(s[0]);

    const seen = new Set<string>();
    return out.filter((c) => c.length <= 3 && !seen.has(c) && (seen.add(c), true));
  };

  return (name: string | null): string | null => {
    if (!name) return null;
    const s = norm(name);
    if (!s) return null;

    const hit = cache.get(s);
    if (hit) return hit;

    const canonical = prot.get(s);
    if (canonical) {
      cache.set(s, canonical);
      return canonical;
    }

    for (const c of candidates(s)) {
      if (!used.has(c)) {
        used.add(c);
        cache.set(s, c);
        return c;
      }
    }
    return null;
  };
};

export const renderChipLabel = (
  rawName: string | null,
  toDisplayName: (name: string | null) => string | null
): React.ReactNode => {
  const raw = rawName || "";
  const displayLabel = isPolyketideMotif(raw) ? raw : toDisplayName(raw) || "X";
  return renderChiralSuperscripts(displayLabel);
};

export const renderTooltipLabel = (
  rawName: string | null,
  toDisplayName: (name: string | null) => string | null
): React.ReactNode => {
  const raw = rawName || "";

  if (isPolyketideMotif(raw)) {
    return toDisplayName(raw) || raw;
  }

  return renderChiralSuperscripts(raw || "Unknown motif");
};