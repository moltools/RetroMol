import { createCookie, getCookie, deleteCookie } from "./utils";

function clearAllCookies() {
  document.cookie.split(";").forEach((c) => {
    const name = c.split("=")[0].trim();
    if (name) document.cookie = `${name}=; expires=Thu, 01 Jan 1970 00:00:00 GMT; path=/;`;
  });
}

describe("cookie utils", () => {
  afterEach(clearAllCookies);

  it("round-trips values containing characters that need URL-encoding", () => {
    // Regression test: getCookie previously returned the raw (still
    // percent-encoded) value instead of decoding it, so anything but a plain
    // alphanumeric value would come back mangled.
    createCookie("greeting", "hello world & friends=you");
    expect(getCookie("greeting")).toBe("hello world & friends=you");
  });

  it("returns null for a cookie that was never set", () => {
    expect(getCookie("nope")).toBeNull();
  });

  it("deleteCookie removes a previously set cookie", () => {
    createCookie("temp", "value");
    expect(getCookie("temp")).toBe("value");

    deleteCookie("temp");
    expect(getCookie("temp")).toBeNull();
  });

  it("cookie remains readable over http (jsdom's default protocol)", () => {
    // Regression test: the Secure attribute must not be set unconditionally,
    // since browsers silently drop Secure cookies set over plain http.
    createCookie("sessionId", "abc-123");
    expect(getCookie("sessionId")).toBe("abc-123");
  });
});
