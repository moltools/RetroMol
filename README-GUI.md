# RetroMol-GUI

Graphical user interface for trying out RetroMol.

This directory contains both a production-ready Docker setup and a developer-friendly local workflow.

## Overview

The system runs five services:
- web: React UI served by nginx
- backend: FLask API served by gunicorn
- db: PostgreSQL with pgvector
- redis: in-memory session and job state store
- maintenance: periodically relabels stale processing jobs

Redis ensures that sessions and job states survive worker restarts and that all backend workers share consistent shared state.

## Build and run with Docker (production mode)

The default setup runs everything containerized:
- Builds and serves the frontend React app behind nginx
- Runs the Flask backend with gunicorn
- Runs an additional backend maintenance script that periodically checks for stale jobs
- Runs PostgreSQL and initializes it from a dump file
- Runs Redis for session/job state
- Exposes a read-only DB user for the backend

### Start the full stack

First make sure to copy `.env.example` to `.env` and adjust any environment variables as needed.

Then run:

```bash
docker compose up -d --build
```

The backend itself loads Redis and DB configuration from `docker/backend.env`.

### Access the application 

- App UI: `http://<server-ip>/**`
- API endpoints: `http://<server-ip>/api/...**`

For local user, `<server-ip>` is typically `localhost:4005`.

### Check container health

Check that the backend and database are reachable through the API:

```bash
curl -i http://<server-ip>/api/health  # should return 200 OK (backend alive)
curl -i http://<server-ip>/api/ready  # should return 200 OK (DB connection OK)
```

For local runs, use:

```bash
curl -i http://localhost:4005/api/health  # should return 200 OK (backend alive)
curl -i http://localhost:4005/api/ready  # should return 200 OK (DB connection OK)
```

> Make sure scripts in `/db/init` are executable before first build:
> ```bash
> chmod +x db/init/*.sh
> ```

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

Then, run the helper script:

```bash
bash ./gui/scripts/dev_backend.sh
```

This script:
- Exports DB_HOST=localhost and REDIS_URL=redis://localhost:6379/0
- Runs Flask in debug mode with auto-reload on port 4000

Verify health endpoint to check backend is running:

```bash
curl -i http://localhost:4000/api/health
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

Local development:
```bash
docker compose -f docker-compose.yml -f docker-compose.dev.yml up -d redis
bash ./gui/scripts/dev_backend.sh
cd ./gui/src/client && npm start
```