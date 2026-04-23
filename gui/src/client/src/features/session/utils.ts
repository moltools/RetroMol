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
};

export function getCookie(name: string): string | null {
  const match = document.cookie
    .split("; ")
    .map(pair => pair.split("="))
    .find(([key]) => key === name);
  return match ? match[1] : null;
};

export function deleteCookie(name: string, path: string = "/"): void {
  document.cookie =
    `${encodeURIComponent(name)}=; ` +
    `expires=Thu, 01 Jan 1970 00:00:00 GMT; ` +
    `path=${path};`;
};