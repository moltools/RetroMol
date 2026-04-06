import { postJson } from "../http";
import { QueryResult, QueryResultSchema, CreateQueryResultRespSchema, GetQueryResultRespSchema } from "./types";
import { z } from "zod";

// Shared option types
export const OrderDirSchema = z.enum(["asc", "desc"]).default("desc");
export type OrderDir = z.infer<typeof OrderDirSchema>;

export const OrderSchema = z.object({
  column: z.string(),  // server whitelists per-query
  dir: OrderDirSchema.optional(),
})
export type Order = z.infer<typeof OrderSchema>;

export const PagingSchema = z.object({
  limit: z.number().int().nonnegative().optional(),
  offset: z.number().int().nonnegative().optional(),
});
export type Paging = z.infer<typeof PagingSchema>;

// Generic payload schema
const QueryPayloadSchema = z.object({
  name: z.string(),
  params: z.record(z.unknown()),
  paging: PagingSchema.optional(),
  order: OrderSchema.optional(),
});
export type QueryPayload = z.infer<typeof QueryPayloadSchema>;

/**
 * Generic query runner. If a `paramSchema` is provided, function will validate
 * the params before sending. Otherwise, it will send params as-is (server will
 * coerce/validate).
 */
export async function runQuery<P extends z.ZodTypeAny = z.ZodRecord<z.ZodString, z.ZodUnknown>>(args: {
  name: string;
  params: z.input<P> | Record<string, unknown>;
  paramSchema?: P;
  paging?: Paging;
  order?: Order;
}): Promise<QueryResult> {
  const safeParams = args.paramSchema
    ? args.paramSchema.parse(args.params)
    : (args.params as Record<string, unknown>);

  const payload: QueryPayload = QueryPayloadSchema.parse({
    name: args.name,
    params: safeParams,
    paging: args.paging,
    order: args.order,
  })

  // Server response is exactly same as QueryResultSchema
  const data = await postJson("/api/query", payload, QueryResultSchema)

  return data;
}