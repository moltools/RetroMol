const pad = (num: number): string => num.toString().padStart(2, "0");

// Formats a duration in seconds as "N day(s) and HH:MM:SS hours".
export function formatUptime(uptimeSeconds: number): string {
  const days = Math.floor(uptimeSeconds / (24 * 3600));
  const hours = Math.floor((uptimeSeconds % (24 * 3600)) / 3600);
  const minutes = Math.floor((uptimeSeconds % 3600) / 60);
  const seconds = Math.floor(uptimeSeconds % 60);
  return `${days} day${days !== 1 ? "s" : ""} and ${pad(hours)}:${pad(minutes)}:${pad(seconds)} hours`;
}
