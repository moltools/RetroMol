import json

from flask import Blueprint, Response, jsonify, request

from routes.session_store import (
    SSE_TICKET_TTL_SECONDS,
    consume_sse_ticket,
    create_sse_ticket,
    load_session_meta,
)

blp_events = Blueprint("events", __name__)
blp_sse_ticket = Blueprint("sse_ticket", __name__)


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
    Compatibility endpoint for browser tabs loaded before session polling.

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

    headers = {
        "Cache-Control": "no-cache",
        "X-Accel-Buffering": "no",
    }
    body = f"event: hello\ndata: {json.dumps({'sessionId': session_id})}\n\n"
    return Response(body, mimetype="text/event-stream", headers=headers)
