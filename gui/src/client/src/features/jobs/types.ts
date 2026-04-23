import type { Session } from "../session/types";
import type { NotificationSeverity } from "../../features/notifications/types";
import { z } from "zod";

export type SetSession = (updater: (prev: Session) => Session) => void;

export const WorkspaceImportDepsSchema = z.object({
  setSession: z.custom<SetSession>(),
  pushNotification: z.function().args(z.string(), z.custom<NotificationSeverity>()).returns(z.void()),
  sessionId: z.string().min(1, "Session ID cannot be empty"),
})

export type WorkspaceImportDeps = z.infer<typeof WorkspaceImportDepsSchema>;

export const BaseNewJobSchema = z.object({
  name: z.string().min(1, "Name cannot be empty"),
})

export const NewCompoundJobSchema = BaseNewJobSchema.extend({
  smiles: z.string().min(1, "SMILES cannot be empty"),
  matchStereochemistry: z.boolean().default(false),
})

export type NewCompoundJob = z.infer<typeof NewCompoundJobSchema>;