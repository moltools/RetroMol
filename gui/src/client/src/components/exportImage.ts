import { toBlob } from "html-to-image";

// Walks up from the captured node to find the first non-transparent background
// -- the node itself is usually a plain, unstyled Box, so without this the
// export would default to a hardcoded color instead of matching the dialog's
// actual (light/dark) paper background.
function resolveBackgroundColor(node: HTMLElement): string {
  let el: HTMLElement | null = node;
  while (el) {
    const bg = getComputedStyle(el).backgroundColor;
    if (bg && bg !== "rgba(0, 0, 0, 0)" && bg !== "transparent") return bg;
    el = el.parentElement;
  }
  return "#ffffff";
}

// Renders `node` (and whatever the user currently has highlighted in it) to a
// downloadable PNG.
export async function exportElementAsPng(node: HTMLElement, filename: string): Promise<void> {
  const blob = await toBlob(node, {
    backgroundColor: resolveBackgroundColor(node),
    cacheBust: true,
    pixelRatio: 2,
    // Fonts are loaded from Google Fonts (see public/index.html); the live DOM
    // already renders with them, so embedding is only needed for fidelity in a
    // *different* environment later. Skipped to avoid a network fetch (and
    // possible failure) on every export.
    skipFonts: true,
    // `node` scrolls horizontally when its content (e.g. a wide backbone
    // structure) overflows the dialog -- capture the full content, not just
    // whatever's currently scrolled into view.
    width: node.scrollWidth,
    height: node.scrollHeight,
    style: { overflow: "visible" },
  });
  if (!blob) throw new Error("Could not render this view to an image.");

  const url = URL.createObjectURL(blob);
  const link = document.createElement("a");
  link.href = url;
  link.download = `${filename}.png`;
  link.click();
  URL.revokeObjectURL(url);
}
