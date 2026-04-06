import { z } from "zod";

export const ItemDrawingResultSchema = z.object({
  svg: z.string(),
});

export type ItemDrawingResult = z.infer<typeof ItemDrawingResultSchema>;