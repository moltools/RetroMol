"""Module for defining session endpoints."""

from flask import Blueprint, request, jsonify

from routes.session_store import (
    MAX_SESSIONS,
    create_session as redis_create_session,
    delete_session as redis_delete_session,
    load_session_with_items,
    merge_session_from_client,
    count_sessions,
    delete_item as redis_delete_item,
    strip_property_from_dict,
)


blp_create_session = Blueprint("create_session", __name__)
blp_delete_session = Blueprint("delete_session", __name__)
blp_get_session = Blueprint("get_session", __name__)
blp_save_session = Blueprint("save_session", __name__)
blp_delete_item = Blueprint("delete_item", __name__)


@blp_create_session.post("/api/createSession")
def create_session() -> tuple[dict[str, str], int]:
    """
    Create a new session.

    :return: a tuple containing a dictionary with the session ID and an HTTP status code.
    """
    payload = request.get_json(force=True) or {}
    new_session = payload.get("session")

    if not isinstance(new_session, dict):
        return {"error": "Missing or invalid session"}, 400

    session_id = new_session.get("sessionId")
    if not isinstance(session_id, str) or not session_id:
        return {"error": "Missing or invalid sessionId"}, 400

    # Hard(ish) cap on total sessions to prevent disk bloat
    current = count_sessions()
    if current >= MAX_SESSIONS:
        return {"error": "Maximum number of sessions reached"}, 503

    # Let Redis throw if something is wrong; we assume new session
    try:
        redis_create_session(new_session)
    except Exception as e:
        return {"error": str(e)}, 400

    return jsonify({"sessionId": session_id}), 200


@blp_delete_session.post("/api/deleteSession")
def delete_session() -> tuple[dict[str, str], int]:
    """
    Delete a session.

    :return: a tuple containing a dictionary with the session ID and an HTTP status code.
    """
    payload = request.get_json(force=True) or {}
    session_id = payload.get("sessionId")

    if not isinstance(session_id, str) or not session_id:
        return {"error": "Missing or invalid sessionId"}, 400

    # Check existence
    sess = load_session_with_items(session_id)
    if sess is None:
        return {"error": "Session not found"}, 404

    redis_delete_session(session_id)

    return jsonify({"sessionId": session_id}), 200


@blp_get_session.post("/api/getSession")
def get_session() -> tuple[dict[str, str], int]:
    """
    Get a session.

    :return: a tuple containing a dictionary with the session data and an HTTP status code.
    """
    payload = request.get_json(force=True) or {}
    session_id = payload.get("sessionId")

    if not isinstance(session_id, str) or not session_id:
        return {"error": "Missing or invalid sessionId"}, 400

    full = load_session_with_items(session_id)
    if full is None:
        return {"error": "Session not found"}, 404

    # Ensure items is always a list so the frontend + Zod don't complain
    items = full.get("items", [])
    if not isinstance(items, list):
        items = []
        full["items"] = items

    # We don't want to send payloads back to the client
    full["items"] = [strip_property_from_dict(item, "payload") for item in items]

    return jsonify({"sessionId": full["sessionId"], "session": full}), 200


@blp_save_session.post("/api/saveSession")
def save_session() -> tuple[dict[str, str], int]:
    """
    Save/merge a session from the client.

    Client owns: which items exist and their user-editable fields
    Server owns: what state those items are in (status, timestamps, errors)

        - for existing items (same id): update fields IN PLACE, but do NOT override
          status/errorMessage/updatedAt (these are owned by the job backend)
        - for new items: accept them as-is
        - items that are missing in the payload are treated as deleted

    :return: a tuple containing a dictionary with the session ID and an HTTP status code.
    """
    payload = request.get_json(force=True) or {}
    new_session = payload.get("session")

    if not isinstance(new_session, dict):
        return {"error": "Missing or invalid session"}, 400

    try:
        merge_session_from_client(new_session)
    except ValueError as e:
        # Use 404 for "Session not found"
        msg = str(e)
        if "Session not found" in msg:
            return {"error": msg}, 404
        return {"error": msg}, 400

    session_id = new_session.get("sessionId")

    return jsonify({"sessionId": session_id}), 200


@blp_delete_item.post("/api/deleteSessionItem")
def delete_item() -> tuple[dict[str, str], int]:
    """
    Delete a single item from a session.

    :return: a tuple containing a dictionary with the operation status and an HTTP status code.
    """
    payload = request.get_json(force=True) or {}
    session_id = payload.get("sessionId")
    item_id = payload.get("itemId")

    if not isinstance(session_id, str) or not session_id:
        return {"error": "Missing or invalid sessionId"}, 400

    if not isinstance(item_id, str) or not item_id:
        return {"error": "Missing or invalid itemId"}, 400

    ok = redis_delete_item(session_id, item_id)
    if not ok:
        return {"error": "Session or item not found"}, 404

    return jsonify({"ok": True}), 200