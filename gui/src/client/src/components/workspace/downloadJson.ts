// Serializes a value to a formatted JSON file and triggers a browser download --
// used for a discovery query item's "Download" button (see WorkspaceDiscovery.tsx)
// so a finished result can be saved and later re-viewed offline via "Load result
// file" (client-side only, see DialogViewDiscoveryQuery.tsx) without ever hitting
// the backend or counting against the saved-query cap.

export function downloadJson(value: unknown, filename: string): void {
  const json = JSON.stringify(value, null, 2);
  const blob = new Blob([json], { type: "application/json;charset=utf-8" });
  const url = URL.createObjectURL(blob);
  const link = document.createElement("a");
  link.href = url;
  link.download = filename;
  document.body.appendChild(link);
  link.click();
  document.body.removeChild(link);
  URL.revokeObjectURL(url);
}
