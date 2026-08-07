"""Shared RQ (Redis Queue) plumbing for every CPU-heavy compute endpoint.

This server has no other background-job mechanism -- every endpoint that actually
computes something (RetroMol parsing, PARAS inference, sequence alignment, RDKit
fingerprinting/conformer search) enqueues its work here rather than running it
in-thread. See routes/jobs.py, routes/discovery.py, and routes/shape.py for the
task functions themselves; this module only holds the queue and the bridge that
lets a synchronous Flask handler use it without any frontend contract change.
"""

import os
import time
from typing import Any, Callable

from prometheus_client import Counter
from redis import Redis
from rq import Queue
from rq.job import Job

REDIS_URL = os.getenv("REDIS_URL", "redis://localhost:6379/0")

# How long a Flask handler blocks waiting for a heavy job before giving up and
# returning 503 -- kept safely under gunicorn's own request timeout (see
# gunicorn.conf.py) so a busy queue degrades to a clear "try again" response
# instead of risking the request itself running long enough to look like a hung
# worker. The RQ job timeout is set a bit longer than this so a job that's still
# legitimately running isn't killed by RQ just because the HTTP caller stopped
# waiting -- it keeps running and its result is simply not observed this time.
HEAVY_JOB_WAIT_TIMEOUT_SECONDS = int(os.getenv("HEAVY_JOB_WAIT_TIMEOUT_SECONDS", "90"))
_HEAVY_JOB_RQ_TIMEOUT_SECONDS = HEAVY_JOB_WAIT_TIMEOUT_SECONDS + 30

_POLL_INTERVAL_SECONDS = 0.25

# A *separate* Redis connection from session_store.redis_client -- RQ wants raw
# bytes, while session_store's connection is configured with decode_responses=True
# for its own JSON-string reads/writes.
_redis = Redis.from_url(REDIS_URL)

queue = Queue("heavy_compute", connection=_redis)

# A separate queue for jobs that are known up front to be slow (currently: a discovery
# query with PMI shape enabled -- see routes/discovery.py's run_discovery_query_job).
# Kept apart from `queue` so a backlog of these can't starve fast jobs (compound
# parsing, cheap queries) -- see docker-compose.yml's `worker`/`worker_light`
# services for how the worker pool is split across the two queues.
queue_pmi = Queue("heavy_compute_pmi", connection=_redis)


def _queue_for(heavy: bool) -> Queue:
    """
    Pick which RQ queue a job should be enqueued on.

    :param heavy: True routes to the protected, PMI/slow-job queue
    :return: the target Queue
    """
    return queue_pmi if heavy else queue

# Labeled by task function name so all migrated endpoints are visible individually
# at /metrics (see app.py's PrometheusMetrics setup) -- exactly the "is the queue
# actually being hit, and how often does it save someone from a slow request"
# question this whole thing started from.
HEAVY_JOB_OUTCOMES = Counter(
    "heavy_job_outcomes_total", "Outcome of heavy_compute job submissions", ["task", "outcome"]
)


class JobStillRunningError(Exception):
    """Raised when a heavy job hasn't finished within the client-facing wait bound."""


def enqueue_job(fn: Callable[..., Any], *args: Any, heavy: bool = False) -> Job:
    """
    Enqueue fn(*args) and return immediately without waiting for it to finish.

    For jobs whose own lifecycle (queued/processing/done/error) is tracked and
    pushed to the client some other way -- e.g. a session item updated via
    routes/session_store.update_item and picked up over SSE -- rather than by the
    HTTP request that submitted it. See enqueue_and_wait for the blocking variant
    most job-submission endpoints use instead.

    :param fn: the task function to run on an RQ worker
    :param args: positional arguments to fn
    :param heavy: True routes this job to the protected, PMI/slow-job queue
    :return: the enqueued RQ Job
    """
    return _queue_for(heavy).enqueue(fn, *args, job_timeout=_HEAVY_JOB_RQ_TIMEOUT_SECONDS)


def enqueue_and_wait(
    fn: Callable[..., Any], *args: Any, timeout: int = HEAVY_JOB_WAIT_TIMEOUT_SECONDS, heavy: bool = False
) -> Any:
    """
    Enqueue fn(*args) and block for up to `timeout` seconds for it to finish,
    returning its result.

    Never blocks past `timeout`: if the job hasn't completed by then (queue
    backlog, all workers busy), raises JobStillRunningError rather than continuing
    to wait -- callers should map that to a 503. `fn` must be a top-level,
    importable function (RQ resolves and pickles it by module path, not by value),
    and both `fn`'s arguments and its return value must be plain, picklable data.

    :param fn: the task function to run on an RQ worker
    :param args: positional arguments to fn
    :param timeout: how long to wait for the result before giving up
    :param heavy: True routes this job to the protected, PMI/slow-job queue
    :return: fn's return value
    :raises JobStillRunningError: if the job hasn't finished within `timeout`
    :raises RuntimeError: if the job itself crashed unexpectedly (not a business-logic
        error -- task functions are expected to catch their own known failure modes
        and return a normal error result, not raise)
    """
    task_name = getattr(fn, "__name__", "unknown")
    job = _queue_for(heavy).enqueue(fn, *args, job_timeout=_HEAVY_JOB_RQ_TIMEOUT_SECONDS)

    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        job.refresh()

        if job.is_finished:
            HEAVY_JOB_OUTCOMES.labels(task=task_name, outcome="completed").inc()
            return job.result

        if job.is_failed:
            HEAVY_JOB_OUTCOMES.labels(task=task_name, outcome="failed").inc()
            raise RuntimeError(f"heavy_compute job failed unexpectedly: {job.exc_info}")

        time.sleep(_POLL_INTERVAL_SECONDS)

    HEAVY_JOB_OUTCOMES.labels(task=task_name, outcome="still_running").inc()
    raise JobStillRunningError(
        "Server is busy processing other requests right now. Please try again in a moment."
    )
