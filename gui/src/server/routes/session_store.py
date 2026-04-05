"""Module for managing session storage using Redis."""

import json
import os
import time
from typing import Any, Callable

import redis

MAX_SESSIONS = 1000
JOB_TIMEOUT_SECONDS = int(os.getenv("JOB_TIMEOUT_SECONDS", "10"))

REDIS_URL = os.getenv("REDIS_URL", "redis://localhost:6379/0")
SESSION_TTL_SECONDS = int(os.getenv("SESSION_TTL_SECONDS", str(7 * 24 * 3600)))

APP_START_KEY = "app:start_epoch"
SESSION_PREFIX = "session:"
ITEM_PREFIX = "session_item:"  # key pattern: session_item:{sessionId}:{itemId}


def _get_redis() -> "redis.Redis":
    """
    Get a Redis client instance.

    :return: Redis client instance
    """
    # Return str isntead of byes (decode_responses=True)
    return redis.Redis.from_url(REDIS_URL, decode_responses=True)


redis_client = _get_redis()


def get_or_init_app_start_epoch() -> int:
    """
    Get or initialize the application start epoch time in Redis.

    :return: the application start epoch time
    """
    now = int(time.time())

    # Try to get existing value; only returns if already set
    if redis_client.setnx(APP_START_KEY, str(now)):
        return now

    # Key exists so read it
    val = redis_client.get(APP_START_KEY)
    try:
        return int(val)
    except ValueError:
        # If corrupted, reset
        redis_client.set(APP_START_KEY, str(now))
        return now


def mark_stale_processing_items() -> int:
    """
    Mark items that have been in 'processing' status for too long as 'error'.

    :return: the number of items updated
    """
    now_ms = int(time.time() * 1000)
    max_age_ms = JOB_TIMEOUT_SECONDS * 1000
    updated_count = 0

    pattern = ITEM_PREFIX + "*"
    for key in redis_client.scan_iter(match=pattern):
        raw = redis_client.get(key)
        if raw is None:
            continue

        try:
            item = json.loads(raw)
        except json.JSONDecodeError:
            continue

        if item.get("status") not in ("processing", "queued"):
            continue

        updated_at = item.get("updatedAt", 0)
        try:
            updated_at = int(updated_at)
        except (TypeError, ValueError):
            updated_at = 0

        age_ms = now_ms - updated_at
        if age_ms > max_age_ms:
            item["status"] = "error"
            item["errorMessage"] = "Processing timed out"
            item["updatedAt"] = now_ms

            redis_client.set(key, json.dumps(item), ex=SESSION_TTL_SECONDS)
            updated_count += 1

    return updated_count


def _session_key(session_id: str) -> str:
    """
    Generate a Redis key for the given session ID.

    :param session_id: the session ID
    :return: the Redis key for the session
    """
    return f"{SESSION_PREFIX}{session_id}"


def _item_key(session_id: str, item_id: str) -> str:
    """
    Generate a Redis key for a session item.

    :param session_id: the session ID
    :param item_id: the item ID
    :return: the Redis key for the session item
    """
    return f"{ITEM_PREFIX}{session_id}:{item_id}"


def count_sessions() -> int:
    """
    Count the number of sessions stored in Redis.

    :return: the number of sessions
    """
    count = 0

    # SESSION_PREFIX already includes "session:"
    pattern = SESSION_PREFIX + "*"
    for _ in redis_client.scan_iter(match=pattern):
        count += 1

    return count


def create_session(session: dict[str, Any]) -> None:
    """
    Create a new session in Redis, storing items separately.

    :param session: the session data to create
    :raises ValueError: if the session does not contain a valid 'sessionId'
    """
    session_id = session.get("sessionId")
    if not isinstance(session_id, str):
        raise ValueError("Session is missing a valid 'sessionId'")

    items = session.get("items", [])
    if not isinstance(items, list):
        items = []

    item_ids: list[str] = []

    # Store each item separately
    for item in items:
        item_id = item.get("id")
        if not item_id:
            continue
        item_ids.append(item_id)
        redis_client.set(
            _item_key(session_id, item_id),
            json.dumps(item),
            ex=SESSION_TTL_SECONDS,
        )

    # Store session metadata; we do NOT embed full items here
    meta = session.copy()
    meta["items"] = item_ids
    redis_client.set(
        _session_key(session_id),
        json.dumps(meta),
        ex=SESSION_TTL_SECONDS,
    )


def load_session_meta(session_id: str) -> dict[str, Any] | None:
    """
    Load session metadata from Redis.

    :param session_id: the session ID to load
    :return: the session metadata, or None if not found
    """
    data = redis_client.get(_session_key(session_id))
    if data is None:
        return None

    return json.loads(data)


def load_session_with_items(session_id: str) -> dict[str, Any] | None:
    """
    Load a full session with all items from Redis.

    :param session_id: the session ID to load
    :return: the full session data, or None if not found
    """
    meta = load_session_meta(session_id)
    if meta is None:
        return None

    item_ids = meta.get("items", []) or []
    if not isinstance(item_ids, list):
        item_ids = []

    items: list[dict] = []
    for item_id in item_ids:
        if not item_id:
            continue
        data = redis_client.get(_item_key(session_id, item_id))
        if data is None:
            continue
        try:
            items.append(json.loads(data))
        except json.JSONDecodeError:
            continue

    # Return full session
    full_session = meta.copy()
    full_session["items"] = items
    return full_session


def delete_session(session_id: str) -> None:
    """
    Delete a session and its items from Redis.

    :param session_id: the session ID to delete
    """
    meta = load_session_meta(session_id)
    if meta is None:
        return

    item_ids = meta.get("items", []) or []
    if not isinstance(item_ids, list):
        item_ids = []

    # Delete each item
    for item_id in item_ids:
        redis_client.delete(_item_key(session_id, item_id))

    # Delete session meta
    redis_client.delete(_session_key(session_id))


def load_item(session_id: str, item_id: str) -> dict[str, Any] | None:
    """
    Load a specific item from a session.

    :param session_id: the session ID
    :param item_id: the item ID
    :return: the item data, or None if not found
    """
    data = redis_client.get(_item_key(session_id, item_id))
    if data is None:
        return None

    return json.loads(data)


def save_item(session_id: str, item: dict[str, Any]) -> None:
    """
    Save a specific item to a session.

    :param session_id: the session ID
    :param item: the item data to save
    """
    item_id = item.get("id")
    if not isinstance(item_id, str):
        raise ValueError("Item is missing a valid 'id'")

    # Ensure the item id is in the session's item list
    meta = load_session_meta(session_id)
    if meta is None:
        raise ValueError(f"Session '{session_id}' does not exist")

    item_ids = meta.get("items", []) or []
    if item_id not in item_ids:
        item_ids.append(item_id)
        meta["items"] = item_ids
        redis_client.set(
            _session_key(session_id),
            json.dumps(meta),
            ex=SESSION_TTL_SECONDS,
        )

    # Save item blob
    redis_client.set(
        _item_key(session_id, item_id),
        json.dumps(item),
        ex=SESSION_TTL_SECONDS,
    )


def update_item(session_id: str, item_id: str, mutator: Callable[[dict[str, Any]], None]) -> bool:
    """
    Update a specific item in a session using a mutator function.

    :param session_id: the session ID
    :param item_id: the item ID
    :param mutator: a function that takes the item dict and modifies it in place
    :return: True if the item was found and updated, False otherwise
    """
    item = load_item(session_id, item_id)
    if item is None:
        return False

    mutator(item)

    save_item(session_id, item)
    return True


# Fields that are owned by the server and should not be overwritten by client data
SERVER_OWNED_FIELDS = {
    "status",
    "errorMessage",
    "retrofingerprints",
    "retrofingerprint512",
    "morganfingerprint2048r2",
    "primarySequences",
    "smiles",
    "taggedSmiles",
    "coverage",
    "updatedAt",
}


def merge_session_from_client(new_session: dict[str, Any]) -> None:
    """
    Merge a session update from the client into the existing session in Redis.

    :param new_session: the session data from the client
    :raises ValueError: if the session does not exist or is invalid
    """
    session_id = new_session.get("sessionId")
    if not isinstance(session_id, str) or not session_id:
        raise ValueError("Session is missing a valid 'sessionId'")

    # Load existing session meta + items
    old_full = load_session_with_items(session_id)
    if old_full is None:
        raise ValueError(f"Session '{session_id}' does not exist")

    old_items_list = old_full.get("items", []) or []
    new_items_list = new_session.get("items", []) or []

    old_ids = {it.get("id") for it in old_items_list if it.get("id")}
    old_by_id: dict[str, dict] = {
        it.get("id"): it
        for it in old_items_list
        if it.get("id")
    }

    merged_items: list[dict] = []
    new_item_ids: list[str] = []

    for new_item in new_items_list:
        item_id = new_item.get("id")
        if not item_id:
            continue

        old_item = old_by_id.get(item_id)

        if old_item is None:
            # New item: accept as-is (client owns everything initially)
            merged_items.append(new_item)
            new_item_ids.append(item_id)
        else:
            # Existing item: merge client fields into old item, presevering server-owned fields
            for key, value in new_item.items():
                if key in SERVER_OWNED_FIELDS:
                    continue
                old_item[key] = value
            merged_items.append(old_item)
            new_item_ids.append(item_id)

    # Delete items that were removed by the client
    new_ids_set = set(new_item_ids)
    deleted_ids = old_ids - new_ids_set
    for item_id in deleted_ids:
        redis_client.delete(_item_key(session_id, item_id))

    # Save merged items
    for item in merged_items:
        save_item(session_id, item)

    # Update session meta, preserving non-item fields
    old_full["items"] = new_item_ids
    # Copy other top-level fields from client session if needed
    for key, value in new_session.items():
        if key == "items":
            continue
        old_full[key] = value

    # Save session meta (without full items array)
    meta = old_full.copy()
    meta["items"] = new_item_ids
    redis_client.set(
        _session_key(session_id),
        json.dumps(meta),
        ex=SESSION_TTL_SECONDS,
    )