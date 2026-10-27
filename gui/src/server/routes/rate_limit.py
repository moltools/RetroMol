"""Per-client rate limiting.

Complements the job queue (routes/queue.py): the queue protects total server
throughput; this protects fairness, so one client can't monopolize it by firing
requests faster than everyone else.
"""

from flask import request
from flask_limiter import Limiter
from prometheus_client import Counter

from routes.queue import REDIS_URL

RATE_LIMIT_REJECTIONS = Counter(
    "rate_limit_rejections_total", "Requests rejected for exceeding a rate limit", ["path"]
)


def _client_key() -> str:
    """
    Identify the calling client for rate-limiting purposes.

    nginx.conf sets `X-Real-IP` to the real client address (`proxy_set_header
    X-Real-IP $remote_addr`) before proxying to this backend, so that's the
    trustworthy value when present. Falls back to `remote_addr` directly for
    dev/direct access that bypasses nginx.

    :return: a string identifying the client
    """
    return request.headers.get("X-Real-IP", request.remote_addr or "unknown")


# storage_uri=REDIS_URL (not an in-process store) so limits are correctly shared
# across every gunicorn worker process and any future backend replica -- an
# in-memory limiter would under-count once more than one process is involved, the
# same correctness gap the old process-local concurrency semaphore had.
limiter = Limiter(
    key_func=_client_key,
    storage_uri=REDIS_URL,
    default_limits=["120 per minute"],
    headers_enabled=True,
)
