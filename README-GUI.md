# RetroMol-GUI

Graphical user interface for trying out RetroMol.

This directory contains both a production-ready Docker setup and a developer-friendly local workflow.

## Overview

The system runs four services:
- **web**: React UI served by nginx (reverse-proxies `/api/*` to backend)
- **backend**: Flask API served by gunicorn (2 worker processes by default, see [Sizing the backend](#sizing-the-backend))
- **worker**: RQ background workers that execute every compute-heavy request (RetroMol parsing, PARAS inference, sequence alignment, RDKit fingerprinting/conformer search) off the `heavy_compute` queue -- 2 replicas by default. Runs the same image as backend, just a different command; see [Background job queue](#background-job-queue)
- **redis**: session/job state store, RQ's job queue broker, and the rate limiter's shared counter store

There is no database *service* -- the compound/BGC database is a read-only DuckDB file, mounted directly into backend and worker (see `RETROMOL_DUCKDB_HOST_PATH` below).

Redis ensures that sessions, job state, and the job queue itself survive individual container restarts, and that backend/worker can scale to multiple replicas while sharing consistent state.

A background maintenance loop (relabeling any session item stuck in `"processing"`, e.g. after a crashed job) runs inside the backend container itself -- started once, in gunicorn's arbiter process, regardless of worker count. There is no separate maintenance container.

## Build and run with Docker (production mode)

The default setup runs everything containerized:
- Builds and serves the frontend React app behind nginx
- Runs the Flask backend with gunicorn
- Runs RQ workers that pick up and execute every heavy-compute request
- Runs Redis for session/job state and the job queue
- Mounts a read-only DuckDB file and the PARAS model file into backend and worker

### Start the full stack

First make sure to copy `.env.example` to `.env` and adjust any environment variables as needed (paths to your DuckDB file and PARAS model, `REDIS_PASSWORD`).

Then run:

```bash
docker compose up -d --build
```

The backend and worker load their Redis/DuckDB/PARAS configuration from `gui/docker/backend.env` (both services share the same environment block in `docker-compose.yml`).

### Access the application

- App UI: `http://<server-ip>/**`
- API endpoints: `http://<server-ip>/api/...**`

For local user, `<server-ip>` is typically `localhost:4005`.

### Check container health

`backend`, `redis`, and `web` each have a Docker healthcheck; `docker compose ps` shows `(healthy)`/`(unhealthy)` directly. `worker` has no HTTP surface to probe -- it relies on `restart: unless-stopped` reacting to the process itself exiting.

You can also check the backend's health endpoints directly:

```bash
curl -i http://<server-ip>/api/health  # 200 OK: the process is up
curl -i http://<server-ip>/api/ready   # 200 OK: DuckDB *and* Redis are both reachable
```

For local runs, use `http://localhost:4005` in place of `<server-ip>`.

### Sizing the backend

Gunicorn worker/thread counts and timeouts live in `gui/src/server/gunicorn.conf.py`, driven by env vars set in `gui/docker/backend.env` -- change them there without rebuilding the image:

| Variable | Default | Meaning |
|---|---|---|
| `GUNICORN_WORKERS` | `2` | gunicorn worker processes |
| `GUNICORN_THREADS` | `4` | threads per worker |
| `GUNICORN_TIMEOUT` | `120` | seconds before gunicorn considers a worker hung |
| `RQ_WORKER_REPLICAS` | `2` | how many `worker` containers process the heavy_compute queue (`docker compose up --scale worker=N` also works if your Compose version doesn't apply `deploy.replicas` outside swarm) |
| `HEAVY_JOB_WAIT_TIMEOUT_SECONDS` | `90` | how long a request blocks waiting on a queued job before returning 503 (kept under `GUNICORN_TIMEOUT`) |

### Background job queue

Every endpoint that actually computes something (compound/gene-cluster parsing, reconstruction, Discovery search/alignment, Tanimoto comparison, PMI shape analysis) enqueues its work on Redis via RQ and blocks briefly for the result, rather than running in the request thread -- the request/response contract is unchanged, but the actual compute happens on the `worker` containers, isolated from the web tier. If the queue is backed up, a request returns `503` with a "try again in a moment" message instead of hanging.

### Observability

- `/metrics` (Prometheus format) is exposed on the backend, including request latency/count by endpoint and custom counters for job outcomes and rate-limit rejections. It is **not** proxied through nginx (only `/api/*` is), so it isn't publicly reachable -- point Prometheus at the `backend` container directly on the Docker network, or `docker exec` in to check it locally:
  ```bash
  docker exec retromol_backend conda run -n retromol-gui --no-capture-output \
    python -c "import urllib.request; print(urllib.request.urlopen('http://localhost:4000/metrics').read().decode())"
  ```
- In production, backend logs are structured JSON (one line per request with method/path/status/duration, plus app events) -- pipe `docker compose logs backend` into `jq` to filter/query them.

### Rate limiting

Per-client request limits (keyed by the `X-Real-IP` header nginx sets) apply on top of a `120/minute` app-wide default: `60/minute` on parsing/reconstruction/Discovery-search endpoints (generous enough for a full batch compound import), `10/minute` on the heavier Tanimoto/PMI-shape comparison endpoints. A breached limit returns `429`.

## Local development mode

You can develop with hot-reloading for both backend and frontend.

### Start only Redis (in Docker)

Expose Redis to your host for local development:

```bash
docker compose -f docker-compose.yml -f docker-compose.dev.yml up -d redis
```

Redis available at:
```
host: localhost
port: 6379
```

### Run the backend locally

First create a virtual environment and install backend dependencies:

```bash
CONDA_SUBDIR=osx-64 conda env create -f ./gui/src/server/environment.backend.dev.yml 
conda activate retromol-gui
```

`environment.backend.dev.yml` only installs the GUI's own dependencies
(`requirements.backend.txt` -- Flask, gunicorn, RQ, etc.). The backend imports
`retromol`, `retromol_alignment`, `retromol_antismash`, `retromol_database`,
`retromol_fingerprint`, and `retromol_synthesis` directly (see e.g.
`gui/src/server/routes/discovery.py`) -- those come from the root package, not from
that env file, so install it in editable mode too, from the repo root (same
`pip install -e /app` step `backend.Dockerfile` runs for the Docker image):

```bash
pip install -e .
```

Then, run the helper script:

```bash
bash ./gui/scripts/dev_backend.sh
```

This script:
- Exports `RETROMOL_DUCKDB_PATH` and `REDIS_URL=redis://localhost:6379/0`
- Runs Flask in debug mode with auto-reload on port 4000

Verify health endpoint to check backend is running:

```bash
curl -i http://localhost:4000/api/health
```

**Also start an RQ worker in a second terminal** (same conda env, same Redis) -- every compute-heavy request (compound/cluster submission, Discovery search, Compare, Shape) now blocks waiting on the `heavy_compute` queue, so without a worker running those requests will just time out after `HEAVY_JOB_WAIT_TIMEOUT_SECONDS` (90s) and return a 503.

The worker runs the task functions itself, so it needs the same environment as the Flask backend (`PYTHONPATH` to import `routes.*`, plus `PARAS_MODEL_PATH`, `CACHE_DIR`, `RETROMOL_DUCKDB_PATH`) -- not just `REDIS_URL`. Use the helper script rather than a bare `rq worker` command:

```bash
conda activate retromol-gui
bash ./gui/scripts/dev_worker.sh
```

By default this one worker listens to both queues (PMI-first, matching a single shared
pile -- see `routes/queue.py`). Optionally, run a second worker terminal dedicated to
the light queue, so a slow PMI-flagged discovery query can never block fast jobs (this
is what production does by default -- see `docker-compose.yml`'s `worker`/`worker_light`
services):

```bash
conda activate retromol-gui
WORKER_QUEUES=heavy_compute bash ./gui/scripts/dev_worker.sh
```

### Run the frontend locally

Make sure to add `.env.development.local` to `src/client` and add the following line for SSE:

```
REACT_APP_SSE_BASE=http://localhost:4000
```

From the React client directory `src/client`, install dependencies and start the development server:

```bash
cd ./gui/src/client
npm install
npm start
```

Ensure the `package.json` has the proxy set to the backend URL:

```json
{
    "proxy": "http://localhost:4000"
}
```

Requests to `/api/...` will automatically proxy to Flask.

## Summary

Production:
```bash
docker compose up -d --build
```

Local development (three terminals: Redis, backend + RQ worker, frontend; run once, before terminal 2: `pip install -e .` from the repo root, in the activated `retromol-gui` conda env):
```bash
docker compose -f docker-compose.yml -f docker-compose.dev.yml up -d redis
bash ./gui/scripts/dev_backend.sh              # terminal 2
bash ./gui/scripts/dev_worker.sh               # terminal 2b, same conda env
cd ./gui/src/client && npm start               # terminal 3
```
