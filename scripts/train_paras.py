#!/usr/bin/env python3
"""
Train (and cache) retromol_paras's PARAS-reimplementation model once, so it can be
mounted read-only into the app instead of every backend/worker container retraining
it from scratch on first use.

Requires hmmpfam2 (HMMER2), hmmscan+hmmpress (HMMER3), and muscle v3.8.1551 on
PATH -- see environment.yml (repo root) for a conda env that installs these
(including an Apple Silicon fallback). Run this with that env active, wherever's
convenient (locally, or via `docker compose run backend ...` since
gui/src/server/environment.backend.yml already has the same tools) -- it doesn't
need to run inside the deployed containers themselves.

Usage:

    python scripts/train_paras.py --cache-dir /path/to/host/paras_cli_cache

Then point docker-compose's PARAS_CACHE_HOST_PATH (see docker-compose.yml's
x-backend-env block) at that same directory, mounted read-only, so the backend/
worker containers find the already-trained model instead of training their own.
Re-run this script (with --force to bypass the cache) and the containers pick up
the new model on their next restart -- nothing else to redeploy.
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

from retromol_paras.train import DEFAULT_CACHE_DIR, train_model


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument(
        "--cache-dir", type=Path, default=DEFAULT_CACHE_DIR,
        help=f"directory to write the trained model + extracted training signatures into (default: {DEFAULT_CACHE_DIR})",
    )
    ap.add_argument("--force", action="store_true", help="retrain from scratch even if a cached model is already present")
    ap.add_argument("-l", "--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"], default="INFO")
    args = ap.parse_args()

    logging.basicConfig(level=args.log_level, format="%(levelname)s %(name)s: %(message)s")

    args.cache_dir.mkdir(parents=True, exist_ok=True)
    model = train_model(cache_dir=args.cache_dir, force=args.force)

    print(f"\nTrained model cached at: {args.cache_dir / 'model.paras.joblib'}")
    print(f"OOB score: {model.oob_score_:.4f}")
    print(f"Substrate classes: {len(model.classes_)}")


if __name__ == "__main__":
    main()
