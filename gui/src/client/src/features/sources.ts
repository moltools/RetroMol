export type EntrySourceLike = {
  name: string;
  databaseName: string;
  url: string | null;
};

export type GroupedSource = {
  databaseName: string;
  count: number;
  items: EntrySourceLike[];
};

/** Group an entry's sources by database, so e.g. two differently-worded MIBiG
 * records for the same compound render as one "MIBiG ×2" chip instead of two
 * identical-looking "MIBiG" chips. */
export function groupSourcesByDatabase(sources: EntrySourceLike[]): GroupedSource[] {
  const byDatabase = new Map<string, EntrySourceLike[]>();
  for (const source of sources) {
    const items = byDatabase.get(source.databaseName) ?? [];
    items.push(source);
    byDatabase.set(source.databaseName, items);
  }
  return Array.from(byDatabase.entries()).map(([databaseName, items]) => ({
    databaseName,
    count: items.length,
    items,
  }));
}
