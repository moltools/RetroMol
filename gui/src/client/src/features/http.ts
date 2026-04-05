import { z } from "zod";

const JSON_HEADERS = { "Content-Type": "application/json" };

export async function postJson<S extends z.ZodTypeAny>(
  url: string,
  body: unknown,
  schema: S
): Promise<z.output<S>> {
  const resp = await fetch(url, {
    method: "POST",
    headers: JSON_HEADERS,
    body: JSON.stringify(body),
  })

  if (!resp.ok) {
    const text = await resp.text().catch(() => "");
    // throw new Error(`${url} failed (${resp.status})${text ? `: ${text}` : ""}`);
    const error: any = new Error(`${url} failed (${resp.status})${text ? `: ${text}` : ""}`);
    error.status = resp.status;
    error.body = text;
    throw error;
  }

  // Try to detect "no body" cases
  const isNoContent =
    resp.status === 204 ||
    resp.headers.get("Content-Length") === "0" ||
    !resp.headers.get("Content-Type")?.includes("application/json");

  const json = isNoContent? undefined : await resp.json().catch(() => {
    throw new Error(`${url} returned non-JSON`)
  })

  const parsed = schema.safeParse(json);
  if (!parsed.success) {
    // Small, actionable error message
    const first = parsed.error.issues[0];
    throw new Error(`${url} response invalid: ${first.path.join(".")}: ${first.message}`);
  }

  return parsed.data;
}