import type { Session } from "../session/types";
import type { NotificationSeverity } from "../notifications/types";
import { z } from "zod";

export const WorkspaceImportDepsSchema = z.object({
  setSession: z.function().args(
    z.function()
      .args(z.custom<Session>())
      .returns(z.custom<Session>())
  )
  .returns(z.void()),
  pushNotification: z.function()
    .args(
      z.string(),
      z.custom<NotificationSeverity>()
    )
    .returns(z.void()),
  sessionId: z.string().min(1, "Session ID cannot be empty"),
})

export type WorkspaceImportDeps = z.infer<typeof WorkspaceImportDepsSchema>;

export const BaseNewJobSchema = z.object({
  name: z.string().min(1, "Name cannot be empty"),
})

export const NewCompoundJobSchema = BaseNewJobSchema.extend({
  smiles: z.string().min(1, "SMILES cannot be empty"),
})

export const NewGeneClusterJobSchema = BaseNewJobSchema.extend({
  fileContent: z.string().min(1, "File content cannot be empty"),
})

export type NewCompoundJob = z.infer<typeof NewCompoundJobSchema>;
export type NewGeneClusterJob = z.infer<typeof NewGeneClusterJobSchema>;