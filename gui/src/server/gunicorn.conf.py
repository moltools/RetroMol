"""Gunicorn configuration.

This is the one place to change worker/thread counts, timeouts, and other server
sizing -- edit the defaults below, or override per-deployment via
GUNICORN_WORKERS/GUNICORN_THREADS/GUNICORN_TIMEOUT in gui/docker/backend.env,
without rebuilding the image.

Also starts the maintenance watchdog loop (see maintenance.py) here, via
on_starting -- that hook runs exactly once in the arbiter process before any
workers are forked, regardless of how many workers are configured, so the loop
runs exactly once no matter what `workers` is set to. Starting it from app.py's
module scope instead would run it once *per worker*, since each worker process
imports app.py independently.
"""

import logging
import os
import threading

bind = f"0.0.0.0:{os.getenv('PORT', '4000')}"
workers = int(os.getenv("GUNICORN_WORKERS", "2"))
threads = int(os.getenv("GUNICORN_THREADS", "4"))
timeout = int(os.getenv("GUNICORN_TIMEOUT", "120"))
accesslog = "-"
errorlog = "-"
loglevel = os.getenv("LOG_LEVEL", "info").lower()


def on_starting(server) -> None:
    """
    Start the maintenance watchdog loop once, in the arbiter process.

    :param server: the gunicorn Arbiter instance (unused, required by gunicorn's hook signature)
    """
    import maintenance

    logging.getLogger(__name__).info(
        "Starting maintenance watchdog thread (interval=%ss)", maintenance.INTERVAL_SECONDS
    )
    threading.Thread(target=maintenance.main, name="maintenance-watchdog", daemon=True).start()
