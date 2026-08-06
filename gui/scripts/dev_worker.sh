#!/usr/bin/env bash
# Run the RQ worker for the heavy_compute queue locally.
# Needs the same environment as dev_backend.sh -- task functions run in this
# process and read PARAS_MODEL_PATH/CACHE_DIR/RETROMOL_DUCKDB_PATH directly,
# and routes.* must be importable here too.
# Usage: ./scripts/dev_worker.sh

set -euo pipefail
cd "$(dirname "$0")/.."  # go to repo root

# --- Setup environment ---
export RETROMOL_DUCKDB_PATH="$HOME/Downloads/retromol.duckdb"

# Redis connection (uses Dockerized Redis)
export REDIS_URL="redis://localhost:6379/0"
export SESSION_TTL_SECONDS=$((7 * 24 * 3600))

# Define cache dir for backend (temp files, etc.)
export CACHE_DIR="$(pwd)/cache"

# Define model paths
export PARAS_MODEL_PATH="$(pwd)/models/all_substrates_model.paras.gz"

# Make sure the worker can import routes.jobs etc.
export PYTHONPATH="$(pwd)/src/server"

echo "Starting RQ worker for queue: heavy_compute"
echo

rq worker --url "${REDIS_URL}" heavy_compute
