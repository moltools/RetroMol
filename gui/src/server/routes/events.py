import json
import time
from typing import Iterator

import redis
from flask import Blueprint, Response, request, stream_with_context

from routes.session_store import REDIS_URL, load_session_meta

blp_events = Blueprint("events", __name__)


def _get_redis() -> "redis.Redis":
    """
    Get a Redis client.

    :return: redis.Redis
    """
    return redis.Redis.from_url(REDIS_URL, decode_responses=True)


@blp_events.get("/api/sessionEvents")
def session_events() -> Response:
    """
    SSE endpoint. Client connects once per session and receives notifications.

    :return: Response
    """
    session_id = request.args.get("sessionId", "")
    if not session_id:
        return Response("missing sessionId\n", status=400, mimetype="text/plain")

    if load_session_meta(session_id) is None:
        return Response("session not found\n", status=404, mimetype="text/plain")

    r = _get_redis()
    pubsub = r.pubsub(ignore_subscribe_messages=True)
    channel = f"session_events:{session_id}"
    pubsub.subscribe(channel)

    def sse(event: str, data: dict) -> str:
        return f"event: {event}\ndata: {json.dumps(data)}\n\n"

    @stream_with_context
    def gen() -> Iterator[str]:
        # Initial hello
        yield sse("hello", {"sessionId": session_id})

        # Keepalive cadence (seconds)
        keepalive_every = 15
        last_keepalive = time.time()

        try:
            while True:
                msg = pubsub.get_message(timeout=1.0)
                if msg and msg.get("type") == "message":
                    try:
                        payload = json.loads(msg["data"])
                    except Exception:
                        payload = {"type": "unknown", "raw": msg.get("data")}

                    # Event name comes from payload["type"]
                    evt_type = payload.get("type", "message")
                    yield sse(evt_type, payload)

                # Keepalive so proxies don't kill the stream
                now = time.time()
                if now - last_keepalive >= keepalive_every:
                    yield sse("keepalive", {"timestamp": int(now * 1000)})
                    last_keepalive = now

        finally:
            try:
                pubsub.close()
            except Exception:
                pass

    headers = {
        "Cache-Control": "no-cache",
        "X-Accel-Buffering": "no",  # helps if behind nginx
        "Connection": "keep-alive",
    }
    return Response(gen(), mimetype="text/event-stream", headers=headers)
