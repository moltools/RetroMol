import { z } from "zod";

export const ServerStartupRespSchema = z.object({
  startup: z.number(),
  current: z.number(),
  uptime: z.number(),
});
export type ServerStartup = z.output<typeof ServerStartupRespSchema>;
