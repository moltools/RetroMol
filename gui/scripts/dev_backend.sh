#!/usr/bin/env bash
# Run the Flask backend locally with auto-reload and local DB connection.
# Usage: ./scripts/dev_backend.sh

set -euo pipefail
cd "$(dirname "$0")/.."  # go to repo root

# --- Setup environment ---
export FLASK_ENV=development
export PORT=4000

# DB connection (uses Dockerized Postgres)
export DB_HOST=localhost
export DB_PORT=5432
export DB_NAME=bionexus
export DB_USER=app_ro
export DB_PASS=apppass_ro

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
echo "DB: ${DB_USER}@${DB_HOST}:${DB_PORT}/${DB_NAME}"
echo

# Run Flask dev server
python -m flask --app app run --host=0.0.0.0 --port="${PORT}" --debug