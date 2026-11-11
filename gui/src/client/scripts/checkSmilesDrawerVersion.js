#!/usr/bin/env node
// Fails if DrawingAttribution.tsx's hardcoded SMILES_DRAWER_VERSION drifts from the
// version npm actually resolved in package-lock.json. Needed because smiles-drawer's
// package.json is blocked from import (its own "exports" field), so the "Drawn with
// SmilesDrawer vX" caption can't read the version at runtime the way ReactionScheme
// does for RDKit -- see DrawingAttribution.tsx for that comparison.
const fs = require("fs");
const path = require("path");

const root = path.join(__dirname, "..");

const lockfile = JSON.parse(fs.readFileSync(path.join(root, "package-lock.json"), "utf8"));
const resolvedVersion = lockfile.packages?.["node_modules/smiles-drawer"]?.version;

if (!resolvedVersion) {
  console.error("checkSmilesDrawerVersion: could not find node_modules/smiles-drawer in package-lock.json");
  process.exit(1);
}

const attributionPath = path.join(root, "src/components/DrawingAttribution.tsx");
const attributionSource = fs.readFileSync(attributionPath, "utf8");
const match = attributionSource.match(/SMILES_DRAWER_VERSION\s*=\s*"([^"]+)"/);

if (!match) {
  console.error(`checkSmilesDrawerVersion: could not find SMILES_DRAWER_VERSION in ${attributionPath}`);
  process.exit(1);
}

const hardcodedVersion = match[1];

if (hardcodedVersion !== resolvedVersion) {
  console.error(
    `checkSmilesDrawerVersion: DrawingAttribution.tsx says smiles-drawer v${hardcodedVersion}, ` +
      `but package-lock.json resolves it to v${resolvedVersion}. ` +
      `Update SMILES_DRAWER_VERSION in src/components/DrawingAttribution.tsx to match.`
  );
  process.exit(1);
}

console.log(`checkSmilesDrawerVersion: OK (v${resolvedVersion})`);
