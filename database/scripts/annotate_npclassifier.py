"""Step 9: classify every compound entry's chemical class via NPClassifier.

Runs after both compound-loading steps (NPAtlas + MIBiG compounds are already
deduplicated by inchikey in `entries` by then -- see load_compounds.py), so each
distinct molecule is classified exactly once no matter which source(s) it came from.
BGC entries are skipped: NPClassifier needs a compound's own SMILES structure
(`entry.raw`), which a bgc entry doesn't have.

Idempotent across reruns via a local JSONL cache (`--cache-path`) keyed by entry_id:
already-classified compounds are skipped without hitting the API again, and only
successes are cached, so a transient failure gets retried on the next run. Requests
are paced (`--requests-per-second`, conservative by default) since NPClassifier is a
free, GNPS2-hosted academic service with no published rate limit.
"""

import argparse
import json
import logging
import time
from pathlib import Path

from npclassifier import ClassificationResult, classify_smiles
from retromol_database.duckdb import RetroMolDuckDB

log = logging.getLogger(__name__)


def _load_cache(cache_path: Path) -> dict[str, ClassificationResult]:
    if not cache_path.exists():
        return {}

    cache: dict[str, ClassificationResult] = {}
    with open(cache_path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            row = json.loads(line)
            cache[row["entry_id"]] = ClassificationResult(
                pathway=row["pathway"],
                superclass=row["superclass"],
                class_=row["class_"],
                is_glycoside=row["is_glycoside"],
            )
    return cache


def _append_cache(cache_path: Path, entry_id: str, result: ClassificationResult) -> None:
    with open(cache_path, "a") as fh:
        fh.write(
            json.dumps(
                {
                    "entry_id": entry_id,
                    "pathway": result.pathway,
                    "superclass": result.superclass,
                    "class_": result.class_,
                    "is_glycoside": result.is_glycoside,
                }
            )
            + "\n"
        )


def _apply(db: RetroMolDuckDB, entry_id: str, result: ClassificationResult) -> None:
    for label in result.pathway:
        db.add_chemical_class_annotation(entry_id, level="pathway", label=label)
    for label in result.superclass:
        db.add_chemical_class_annotation(entry_id, level="superclass", label=label)
    for label in result.class_:
        db.add_chemical_class_annotation(entry_id, level="class", label=label)
    if result.is_glycoside:
        db.add_chemical_class_annotation(entry_id, level="is_glycoside", label="Yes")


def run(
    db_path: str | Path,
    cache_path: str | Path,
    requests_per_second: float = 2.0,
    limit: int | None = None,
    log_every: int = 100,
) -> None:
    cache_path = Path(cache_path)
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache = _load_cache(cache_path)
    log.info("annotate_npclassifier: loaded %d cached classifications from %s", len(cache), cache_path)

    min_interval = 1.0 / requests_per_second if requests_per_second > 0 else 0.0
    last_call = 0.0

    classified = 0
    reused = 0
    skipped_no_smiles = 0
    failed = 0

    db = RetroMolDuckDB.open(db_path)
    try:
        for entry in db.iter_entries():
            if entry.type != "compound":
                continue

            if entry.id in cache:
                _apply(db, entry.id, cache[entry.id])
                reused += 1
                continue

            if not entry.raw:
                skipped_no_smiles += 1
                continue

            if limit is not None and classified >= limit:
                continue

            elapsed = time.monotonic() - last_call
            if elapsed < min_interval:
                time.sleep(min_interval - elapsed)
            last_call = time.monotonic()

            result = classify_smiles(entry.raw)
            if result is None:
                failed += 1
                continue

            _apply(db, entry.id, result)
            _append_cache(cache_path, entry.id, result)
            classified += 1

            if log_every > 0 and classified % log_every == 0:
                log.info(
                    "annotate_npclassifier: classified=%d reused=%d failed=%d skipped_no_smiles=%d",
                    classified, reused, failed, skipped_no_smiles,
                )
    finally:
        db.close()

    log.info(
        "annotate_npclassifier: classified=%d reused=%d failed=%d skipped_no_smiles=%d",
        classified, reused, failed, skipped_no_smiles,
    )


def main() -> None:
    logging.basicConfig(level=logging.INFO)

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--db-path", required=True)
    ap.add_argument("--cache-path", required=True)
    ap.add_argument("--requests-per-second", type=float, default=2.0)
    ap.add_argument("--limit", type=int, default=None, help="classify at most N new compounds (testing/dry-run)")
    ap.add_argument("--log-every", type=int, default=100)
    args = ap.parse_args()

    run(
        db_path=args.db_path,
        cache_path=args.cache_path,
        requests_per_second=args.requests_per_second,
        limit=args.limit,
        log_every=args.log_every,
    )


if __name__ == "__main__":
    main()
