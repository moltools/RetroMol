#!/usr/bin/env bash
# Run the Flask backend locally with auto-reload and local DB connection.
# Usage: ./scripts/dev_backend.sh

set -euo pipefail
cd "$(dirname "$0")/.."  # go to repo root

# --- Setup environment ---
export FLASK_ENV=development
export PORT=4000

# DB connection
export RETROMOL_DUCKDB_PATH="$HOME/Downloads/retromol.duckdb"

# Redis connection (uses Dockerized Redis)
export REDIS_URL="redis://localhost:6379/0"
export SESSION_TTL_SECONDS=$((7 * 24 * 3600))

# Define cache dir for backend (temp files, etc.)
export CACHE_DIR="$(pwd)/cache"

# Define model paths
export PARAS_MODEL_PATH="$(pwd)/models/all_substrates_model.paras.gz"

# Make sure Flask can find the app
export PYTHONPATH="$(pwd)/src/server"

echo "Starting Flask backend on http://localhost:${PORT} (hot reload enabled)"
echo "DuckDB: ${RETROMOL_DUCKDB_PATH}"
echo

# Run Flask dev server
python -m flask --app app run --host=0.0.0.0 --port="${PORT}" --debug