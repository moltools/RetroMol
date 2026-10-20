import json
import time
from typing import Iterator

import redis
from flask import Blueprint, Response, jsonify, request, stream_with_context

from routes.session_store import (
    REDIS_URL,
    SSE_TICKET_TTL_SECONDS,
    consume_sse_ticket,
    create_sse_ticket,
    load_session_meta,
)

blp_events = Blueprint("events", __name__)
blp_sse_ticket = Blueprint("sse_ticket", __name__)


def _get_redis() -> "redis.Redis":
    """
    Get a Redis client.

    :return: redis.Redis
    """
    return redis.Redis.from_url(REDIS_URL, decode_responses=True)


@blp_sse_ticket.post("/api/sessionEventsTicket")
def session_events_ticket() -> tuple[Response, int]:
    """
    Mint a short-lived, single-use ticket for opening the SSE stream.

    The frontend calls this (a normal POST, so sessionId stays out of the
    URL/logs) and passes the returned ticket as the ?ticket= query param when
    opening the EventSource connection below.

    :return: a tuple containing a dictionary with the ticket and an HTTP status code
    """
    payload = request.get_json(force=True) or {}
    session_id = payload.get("sessionId")

    if not isinstance(session_id, str) or not session_id:
        return jsonify({"error": "Missing or invalid sessionId"}), 400

    if load_session_meta(session_id) is None:
        return jsonify({"error": "Session not found"}), 404

    ticket = create_sse_ticket(session_id)
    return jsonify({"ticket": ticket, "expiresInSeconds": SSE_TICKET_TTL_SECONDS}), 200


@blp_events.get("/api/sessionEvents")
def session_events() -> Response:
    """
    SSE endpoint. Client connects once per session and receives notifications.

    Authorized via a one-time ticket (see /api/sessionEventsTicket above)
    rather than a raw sessionId, since EventSource can only send a GET and
    the sessionId must not end up sitting in a URL.

    :return: Response
    """
    ticket = request.args.get("ticket", "")
    if not ticket:
        return Response("missing ticket\n", status=400, mimetype="text/plain")

    session_id = consume_sse_ticket(ticket)
    if session_id is None:
        return Response("invalid or expired ticket\n", status=401, mimetype="text/plain")

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
