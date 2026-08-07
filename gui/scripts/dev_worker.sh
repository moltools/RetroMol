#!/usr/bin/env bash
# Run an RQ worker locally. Needs the same environment as dev_backend.sh -- task
# functions run in this process and read PARAS_MODEL_PATH/CACHE_DIR/
# RETROMOL_DUCKDB_PATH directly, and routes.* must be importable here too.
#
# By default listens to both queues (heavy_compute_pmi first, then heavy_compute),
# matching a single local worker sharing one pile across everything -- see
# routes/queue.py. To test the production split locally (a worker dedicated to
# light jobs that never blocks behind a PMI job), run a second terminal with:
#   WORKER_QUEUES=heavy_compute bash gui/scripts/dev_worker.sh
#
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

QUEUES="${WORKER_QUEUES:-heavy_compute_pmi heavy_compute}"

echo "Starting RQ worker for queue(s): ${QUEUES}"
echo

rq worker --url "${REDIS_URL}" ${QUEUES}
