import React from "react";
import { render } from "@testing-library/react";
import { SvgViewer } from "./SvgViewer";

describe("SvgViewer", () => {
  it("strips <script> tags and inline event handlers from untrusted SVG before rendering", () => {
    // Regression test for the dangerouslySetInnerHTML XSS finding: the
    // backend returns SVG derived from user-uploaded compounds/gene
    // clusters, so it must never be injected verbatim.
    const maliciousSvg =
      '<svg onload="window.__pwned = true"><script>window.__pwned = true;</script><circle cx="5" cy="5" r="4" /></svg>';

    const { container } = render(<SvgViewer svg={maliciousSvg} />);

    // This asserts the *absence* of dangerous DOM nodes/attributes, which
    // has no accessible-role equivalent to query by — direct node access is
    // the correct tool here, not a shortcut around Testing Library.
    /* eslint-disable testing-library/no-container, testing-library/no-node-access */
    expect(container.querySelector("script")).toBeNull();
    expect(container.querySelector("svg")?.hasAttribute("onload")).toBe(false);
    // legitimate markup still survives sanitization
    expect(container.querySelector("circle")).not.toBeNull();
    /* eslint-enable testing-library/no-container, testing-library/no-node-access */
  });
});
