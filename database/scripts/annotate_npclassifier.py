"""Step 9: classify every compound entry's chemical class via NPClassifier.

Runs after both compound-loading steps (NPAtlas + MIBiG compounds are already
deduplicated by inchikey in `entries` by then -- see load_compounds.py), so each
distinct molecule is classified exactly once no matter which source(s) it came from.
BGC entries are skipped: NPClassifier needs a compound's own SMILES structure
(`entry.raw`), which a bgc entry doesn't have.

Idempotent across reruns via a local JSONL cache (`--cache-path`) keyed by entry_id:
already-classified compounds are skipped without hitting the API again, and only
successes are cached -- flushed to disk as each one completes, not batched at the end
-- so a kill mid-run loses at most the handful of requests in flight, and a transient
failure gets retried on the next run.

Classification is I/O-bound (waiting on NPClassifier's API), not CPU-bound, so
`--workers` runs multiple requests concurrently via a thread pool -- the GIL doesn't
block this the way it would CPU-bound work, since each thread is blocked on network
I/O rather than holding the GIL. `--requests-per-second` still caps the *combined*
rate across all workers (a shared, thread-safe RateLimiter -- see below), since that
cap is about being a good citizen towards NPClassifier's free, GNPS2-hosted service,
not about this process's own resources. DB writes and cache appends only ever happen
on the main thread (as results complete), so there's no concurrent-write concern
there -- only the HTTP calls themselves run in parallel.
"""

import argparse
import itertools
import json
import logging
import threading
import time
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from pathlib import Path

from tqdm import tqdm

from npclassifier import ClassificationResult, classify_smiles
from retromol_database.duckdb import Entry, RetroMolDuckDB

log = logging.getLogger(__name__)


class RateLimiter:
    """Thread-safe request pacing shared across worker threads.

    Each call reserves the next available time slot under a lock (cheap, no waiting
    while holding it), then sleeps outside the lock until its own slot arrives -- so
    threads don't serialize on the actual waiting, only on the quick slot reservation.
    """

    def __init__(self, requests_per_second: float) -> None:
        self._min_interval = 1.0 / requests_per_second if requests_per_second > 0 else 0.0
        self._lock = threading.Lock()
        self._next_slot = time.monotonic()

    def wait_for_slot(self) -> None:
        with self._lock:
            now = time.monotonic()
            slot = max(now, self._next_slot)
            self._next_slot = slot + self._min_interval
        delay = slot - now
        if delay > 0:
            time.sleep(delay)


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
    workers: int = 8,
    limit: int | None = None,
    log_every: int = 100,
) -> None:
    cache_path = Path(cache_path)
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache = _load_cache(cache_path)
    log.info("annotate_npclassifier: loaded %d cached classifications from %s", len(cache), cache_path)

    rate_limiter = RateLimiter(requests_per_second)

    classified = 0
    reused = 0
    skipped_no_smiles = 0
    failed = 0

    def classify_one(entry: Entry) -> tuple[Entry, ClassificationResult | None]:
        rate_limiter.wait_for_slot()
        return entry, classify_smiles(entry.raw)

    db = RetroMolDuckDB.open(db_path)
    try:
        to_classify: list[Entry] = []
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

            if limit is not None and len(to_classify) >= limit:
                continue

            to_classify.append(entry)

        log.info(
            "annotate_npclassifier: %d compounds to classify (workers=%d, requests_per_second=%s)",
            len(to_classify), workers, requests_per_second,
        )

        # Bounded sliding window rather than submitting the whole backlog up front --
        # at hundreds/thousands of compounds, an upfront submit() for everything means
        # a Ctrl+C has to wait for every already-launched request (each with its own
        # retry/backoff chain, worse under rate-limiting) before the pool can actually
        # exit, since ThreadPoolExecutor won't drop already-running work. Keeping at
        # most `max_pending` in flight means an interrupt only has to wait for that many
        # -- same shape as parse_gbks.py's own bounded-window loop, same reason.
        max_pending = max(workers * 2, 1)
        entries_iter = iter(to_classify)

        pool = ThreadPoolExecutor(max_workers=workers)
        try:
            pending = {pool.submit(classify_one, e) for e in itertools.islice(entries_iter, max_pending)}

            with tqdm(total=len(to_classify), desc="annotate_npclassifier", unit="cmpd") as pbar:
                while pending:
                    done_futures, pending = wait(pending, return_when=FIRST_COMPLETED)

                    for future in done_futures:
                        entry, result = future.result()

                        if result is None:
                            failed += 1
                        else:
                            _apply(db, entry.id, result)
                            _append_cache(cache_path, entry.id, result)
                            classified += 1

                        pbar.update(1)
                        pbar.set_postfix(classified=classified, reused=reused, failed=failed)

                        done = classified + failed
                        if log_every > 0 and done % log_every == 0:
                            log.info(
                                "annotate_npclassifier: classified=%d reused=%d failed=%d skipped_no_smiles=%d",
                                classified, reused, failed, skipped_no_smiles,
                            )

                        next_entry = next(entries_iter, None)
                        if next_entry is not None:
                            pending.add(pool.submit(classify_one, next_entry))
        except KeyboardInterrupt:
            log.warning("annotate_npclassifier: interrupted -- cancelling not-yet-started requests")
            pool.shutdown(wait=False, cancel_futures=True)
            raise
        else:
            pool.shutdown(wait=True)
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
    ap.add_argument("--workers", type=int, default=8, help="concurrent requests (I/O-bound, not CPU-bound)")
    ap.add_argument("--limit", type=int, default=None, help="classify at most N new compounds (testing/dry-run)")
    ap.add_argument("--log-every", type=int, default=100)
    args = ap.parse_args()

    run(
        db_path=args.db_path,
        cache_path=args.cache_path,
        requests_per_second=args.requests_per_second,
        workers=args.workers,
        limit=args.limit,
        log_every=args.log_every,
    )


if __name__ == "__main__":
    main()
