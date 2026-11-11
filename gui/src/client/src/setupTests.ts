// jest-dom adds custom jest matchers for asserting on DOM nodes.
// https://github.com/testing-library/jest-dom
import "@testing-library/jest-dom";

// The jsdom version bundled with this react-scripts/jest setup doesn't
// implement the Web Crypto API, but app code calls crypto.randomUUID()
// directly (it's available in every real browser this app targets).
if (typeof globalThis.crypto === "undefined" || typeof globalThis.crypto.randomUUID !== "function") {
  // eslint-disable-next-line @typescript-eslint/no-var-requires
  const { webcrypto } = require("node:crypto");
  Object.defineProperty(globalThis, "crypto", { value: webcrypto });
}
