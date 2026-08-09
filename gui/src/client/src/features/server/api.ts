import { getJson } from "../http";
import { ServerStartup, ServerStartupRespSchema } from "./types";

export async function getServerStartup(signal?: AbortSignal): Promise<ServerStartup> {
  return getJson("/api/startup", ServerStartupRespSchema, signal);
}
