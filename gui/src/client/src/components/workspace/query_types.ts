export type SequenceItem = {
  id: string;
  isGap: boolean;
  name: string | null;
  smiles: string | null;
};

export type Sequence = {
  id: string;
  name: string | null;
  sequence: SequenceItem[];
};

export type MsaRow = {
  id: string;
  name?: string;
  kind?: "compound" | "cluster" | null;
  db_id?: number | null;
  alignment_score: number | null;
  cosine_score: number | null;
  match_score: number | null;
  sequence: Sequence[];
};

export type QueryResult = {
  msa: MsaRow[];
};