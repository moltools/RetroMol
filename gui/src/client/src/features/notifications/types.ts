import { z } from "zod";

export const NotificationSeveritySchema = z.enum([
  "info",
  "warning",
  "error",
  "success",
]);

export type NotificationSeverity = z.infer<typeof NotificationSeveritySchema>;