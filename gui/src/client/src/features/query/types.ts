import { z } from "zod";

export const QueryResultSchema = z.object({
  name: z.string(),
  columns: z.array(z.string()),
  rows: z.array(z.record(z.unknown())),
  rowCount: z.number().int().nonnegative(),
  limit: z.number().int().nonnegative(),
  offset: z.number().int().nonnegative(),
  elapsed_ms: z.number().nonnegative(),
});

export type QueryResult = z.infer<typeof QueryResultSchema>;

// Simple response wrapper
export const CreateQueryResultRespSchema = z.object({ result: QueryResultSchema });
export const GetQueryResultRespSchema = z.object({ result: QueryResultSchema });