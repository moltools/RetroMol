import React from "react";
import Typography from "@mui/material/Typography";
import Chip from "@mui/material/Chip";
import Alert from "@mui/material/Alert";
import CircularProgress from "@mui/material/CircularProgress";
import { DialogWindow } from "../DialogWindow";
import { MsaRow } from "./query_types";

type DialogRowInfoProps = {
  open: boolean;
  onClose: () => void;
  msaRow: MsaRow;
};

type Annotation = {
  scheme: string;
  key: string;
  value: string;
};

type Reference = {
  name: string;
  database_name: string;
  database_identifier: string;
};

type RowInfoResponse = {
  annotations: Annotation[];
  references: Reference[];
};

function referenceToUrl(ref: Reference): string | null {
  const { database_name, database_identifier } = ref;

  switch (database_name.toLowerCase()) {
    case "mibig":
      return `https://mibig.secondarymetabolites.org/repository/${database_identifier}`;
    case "npatlas":
      return `https://www.npatlas.org/explore/compounds/${database_identifier}`
    default:
      return null;
  };
};

function referenceLabel(ref: Reference): string {
  return `${ref.name} (${ref.database_name}: ${ref.database_identifier})`;
};

export const DialogRowInfo: React.FC<DialogRowInfoProps> = ({
  open,
  onClose,
  msaRow,
}) => {
  const kind = msaRow.kind ?? null;
  const dbId = msaRow.db_id ?? null;

  const canFetch = open && kind != null && dbId != null;

  const [loading, setLoading] = React.useState(false);
  const [error, setError] = React.useState<string | null>(null);
  const [data, setData] = React.useState<RowInfoResponse | null>(null);

  const cacheRef = React.useRef(new Map<string, RowInfoResponse>());

  const cacheKey = kind && dbId != null ? `${kind}:${dbId}` : null;

  const sortedAnnotations = React.useMemo(() => {
    if (!data?.annotations) return [];
    return [...data.annotations].sort((a, b) => {
      const s = a.scheme.localeCompare(b.scheme);
      if (s !== 0) return s;

      const k = a.key.localeCompare(b.key);
      if (k !== 0) return k;

      return a.value.localeCompare(b.value);
    });
  }, [data?.annotations]);

  const sortedReferences = React.useMemo(() => {
    if (!data?.references) return [];
    return [...data.references].sort((a, b) => {
      const d = a.database_name.localeCompare(b.database_name);
      if (d !== 0) return d;

      const i = a.database_identifier.localeCompare(b.database_identifier);
      if (i !== 0) return i;

      return a.name.localeCompare(b.name);
    });
  }, [data?.references]);

  React.useEffect(() => {
    if (!canFetch || !cacheKey) {
      setLoading(false);
      setError(null);
      setData(null);
      return;
    }

    // Serve from cache immediately if present
    const cached = cacheRef.current.get(cacheKey);
    if (cached) {
      setLoading(false);
      setError(null);
      setData(cached);
      return;
    }

    // Otherwise fetch
    const controller = new AbortController();
    let alive = true;

    (async () => {
      try {
        setLoading(true);
        setError(null);

        const res = await fetch(
          `/api/itemInfo?kind=${encodeURIComponent(kind!)}&db_id=${dbId}`,
          {
            method: "GET",
            signal: controller.signal,
            headers: { Accept: "application/json" },
          }
        );

        if (!res.ok) {
          throw new Error(`Error fetching data: ${res.status} ${res.statusText}`);
        }

        const json = (await res.json()) as RowInfoResponse;
        if (!alive) return;

        // Store in cache
        cacheRef.current.set(cacheKey, json);

        setData(json);
      } catch (e: any) {
        if (e?.name === "AbortError") return;
        if (!alive) return;
        setError(e?.message ?? "Unknown error");
      } finally {
        if (!alive) return;
        setLoading(false);
      }
    })();

    return () => {
      alive = false;
      controller.abort();
    };
  }, [canFetch, cacheKey, kind, dbId]);

  return (
    <DialogWindow
      open={open}
      onClose={onClose}
      title="Row information"
      dividers
      actions={[
        { label: "Close", variant: "text", color: "inherit", onClick: onClose },
      ]}
    >
      {(!canFetch || (!data?.annotations.length && !data?.references.length)) ? (
        <Typography variant="body2">
          No additional information available for this row.
        </Typography>
      ) : null}

      {data?.references.length ? (
        <>
          <Typography variant="h6" gutterBottom>
            References
          </Typography>
          {sortedReferences.map((reference) => {
            const url = referenceToUrl(reference);

            return (
              <Chip
                key={`${reference.database_name}:${reference.database_identifier}:${reference.name}`}
                label={referenceLabel(reference)}
                clickable={Boolean(url)}
                component={url ? "a" : "div"}
                href={url ?? undefined}
                target={url ? "_blank" : undefined}
                rel={url ? "noopener noreferrer" : undefined}
                sx={{
                  m: 0.5,
                  ...(url && {
                    border: "1px solid",
                    borderColor: "primary.main",
                  }),
                }}
              />
            );
          })}
        </>
      ) : null}

      {data?.annotations.length ? (
        <>
          <Typography variant="h6" gutterBottom style={{ marginTop: "16px" }}>
            Annotations
          </Typography>
          {sortedAnnotations.map((annotation, index) => (
            <Chip
              key={index}
              label={`${annotation.scheme}: ${annotation.key} = ${annotation.value}`}
              style={{ margin: "4px" }}
            />
          ))}
        </>
      ) : null}

      {loading && (
        <CircularProgress size={24} />
      )}

      {error && (
        <Alert severity="error">
          {error}
        </Alert>
      )}
    </DialogWindow>
  );
};