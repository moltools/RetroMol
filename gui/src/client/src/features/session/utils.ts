export function createCookie(name: string, value: string, days?: number, path: string = "/"): void {
  let cookieStr = `${encodeURIComponent(name)}=${encodeURIComponent(value)}`;
  if (typeof days === "number") {
    const expires = new Date();
    // Set the expiration date for the cookie; here we set it to "days" days from now
    expires.setTime(expires.getTime() + days * 24 * 60 * 60 * 1000);
    cookieStr += `; expires=${expires.toUTCString()}`;
  }
  cookieStr += `; path=${path}; samesite=strict`;
  // "Secure" cookies are silently dropped by the browser over plain http, so
  // only add it when we're actually on https (keeps local http dev working).
  if (window.location.protocol === "https:") {
    cookieStr += "; secure";
  }
  document.cookie = cookieStr;
};

export function getCookie(name: string): string | null {
  const match = document.cookie
    .split("; ")
    .map(pair => pair.split("="))
    .find(([key]) => key === name);
  return match ? decodeURIComponent(match[1]) : null;
};

export function deleteCookie(name: string, path: string = "/"): void {
  document.cookie =
    `${encodeURIComponent(name)}=; ` +
    `expires=Thu, 01 Jan 1970 00:00:00 GMT; ` +
    `path=${path};`;
};