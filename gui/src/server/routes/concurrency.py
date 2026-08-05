"""Process-wide concurrency guard for CPU-heavy compute endpoints.

This server has no background job queue -- CPU-heavy work (conformer search,
RDKit fingerprinting) runs synchronously inside one of gunicorn's shared request
threads (see gui/docker/backend.Dockerfile: -w 1 --threads 4, and routes/jobs.py's
compound-parsing flow, which follows the same pattern). Every endpoint in the app
shares that same small thread pool. Without a cap, enough concurrent heavy requests
can saturate every thread at once; if the single worker process then goes quiet long
enough, gunicorn's --timeout kills and restarts it -- a full outage for every
connected user, not just the ones running heavy computations.

This guard fails fast rather than letting requests pile up toward that failure mode:
once the cap is reached, further requests get an immediate 503 asking the client to
retry, instead of blocking/queueing for a thread that may not free up in time.

Process-local only -- coordinates concurrency within a single backend replica. If
this backend is ever scaled to more than one replica, this needs to move to a
Redis-backed counter (or a real job queue) to stay globally accurate.
"""

import os
import threading
from contextlib import contextmanager
from typing import Iterator

# Deliberately well under the thread pool size (4) so heavy computations can never
# occupy every thread -- lightweight endpoints (session saves, discovery lookups,
# health checks) always have room to run.
HEAVY_COMPUTE_MAX_CONCURRENCY = int(os.getenv("HEAVY_COMPUTE_MAX_CONCURRENCY", "2"))

_slots = threading.BoundedSemaphore(HEAVY_COMPUTE_MAX_CONCURRENCY)


class ServerBusyError(Exception):
    """Raised when the heavy-compute concurrency cap is already saturated."""


@contextmanager
def heavy_compute_slot() -> Iterator[None]:
    """
    Reserve one of the limited heavy-compute slots for the duration of the block.

    Never blocks: if none are free, raises ServerBusyError immediately rather than
    waiting for one to open up.

    :raises ServerBusyError: if all HEAVY_COMPUTE_MAX_CONCURRENCY slots are taken
    """
    if not _slots.acquire(blocking=False):
        raise ServerBusyError(
            f"Server is busy running {HEAVY_COMPUTE_MAX_CONCURRENCY} other shape/similarity "
            "computation(s) right now. Please try again in a moment."
        )
    try:
        yield
    finally:
        _slots.release()
