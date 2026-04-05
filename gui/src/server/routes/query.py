"""Module for defining database query endpoints."""

import os
import time

import psycopg
from flask import Blueprint, request, jsonify
from psycopg import sql
from pgvector.psycopg import register_vector

from routes.query_registry import QUERIES

blp = Blueprint("query", __name__)

# Knobs
DEFAULT_LIMIT = 500
MAX_OFFSET = 50_000
STATEMENT_TIMEOUT_MS = 3000  # 3 seconds


def dsn_from_env() -> str:
    """
    Construct the Postgres DSN from environment variables.

    :return: the Postgres DSN string
    """
    dsn = os.getenv("DATABASE_URL")
    if dsn:
        return dsn
    host = os.getenv("DB_HOST", "db")
    port = os.getenv("DB_PORT", "5432")
    name = os.getenv("DB_NAME", "bionexus")
    user = os.getenv("DB_USER", "app_ro")
    pwd = os.getenv("DB_PASS") or os.getenv("DB_PASSWORD", "apppass_ro")
    return f"postgresql://{user}:{pwd}@{host}:{port}/{name}"


def coerce_params(spec: dict[str, type], data: dict) -> tuple[dict, str | None]:
    """
    Coerce and validate parameters from the request data according to the spec.

    :param spec: a dictionary mapping parameter names to expected types
    :param data: the input data dictionary
    :return: a tuple containing the coerced parameters dictionary and an error message (or None
    """
    out = {}
    for k, typ in spec.items():
        if k not in data:
            return {}, f"Missing param: {k}"
        v = data[k]
        try:
            if typ is float:
                out[k] = float(v)
            elif typ is int:
                out[k] = int(v)
            elif typ is str:
                out[k] = str(v)
            else:
                out[k] = v
        except Exception:
            return {}, f"Invalid type for {k}"
    return out, None


def execute_named_query(
        name: str,
        params: dict | None = None,
        paging: dict | None = None,
        order: dict | None = None,
) -> dict:
    """
    Execute a predefined database query with parameters, paging, and ordering.

    :param name: the name of the predefined query
    :param params: a dictionary of query parameters
    :param paging: a dictionary with 'limit' and 'offset' for paging
    :param order: a dictionary with 'column' and 'dir' for ordering
    :return: a dictionary with query results
    :raises ValueError: if there is a parameter validation error
    :raises TimeoutError: if the query times out
    :raises RuntimeError: if there is a database error
    """
    params = params or {}
    paging = paging or {}
    order = order or {}

    if not name or name not in QUERIES:
        raise ValueError("Invalid or missing query name")

    qinfo = QUERIES[name]

    # Validate required/optional params
    req_spec = qinfo.get("required", {})
    opt_spec = qinfo.get("optional", {})
    typed, err = coerce_params(req_spec, params)
    if err:
        raise ValueError(err)
    # Optional params (coerce if present)
    for k, typ in opt_spec.items():
        if k in params:
            try:
                if typ is float:
                    typed[k] = float(params[k])
                elif typ is int:
                    typed[k] = int(params[k])
                elif typ is str:
                    typed[k] = str(params[k])
                else:
                    typed[k] = params[k]
            except Exception:
                raise ValueError(f"Invalid type for {k}")

    # Preprocess params for specific queries
    preprocess = qinfo.get("preprocess_params")
    if preprocess:
        try:
            typed = preprocess(typed) or typed
        except Exception as e:
            raise ValueError(f"Parameter preprocessing error: {str(e)}")

    # Paging
    limit = int(paging.get("limit", DEFAULT_LIMIT))
    offset = int(paging.get("offset", 0))
    offset = max(0, min(MAX_OFFSET, offset))
    typed["limit"] = limit
    typed["offset"] = offset

    # Order-by (whitelisted)
    allowed_cols = qinfo.get("allowed_order_cols", set())
    order_col = (order.get("column") if isinstance(order, dict) else None) or qinfo.get("default_order_col")
    if order_col not in allowed_cols:
        order_col = qinfo.get("default_order_col")
    order_dir = (order.get("dir", qinfo.get("default_order_dir", "ASC")) if isinstance(order, dict) else qinfo.get(
        "default_order_dir", "ASC"))
    order_dir = "DESC" if str(order_dir).upper().startswith("D") else "ASC"

    # Render final SQL safely for the ORDER BY identifier
    base_sql = qinfo["sql"]
    rendered = (
            base_sql.format(
                order_col=sql.Identifier(order_col).as_string(psycopg.connect(dsn_from_env())),
                order_dir=order_dir,
            ).rstrip().rstrip(";")
            + f" LIMIT %(limit)s OFFSET %(offset)s"
    )

    # Exec (read-only, short timeout, public schema)
    dsn = dsn_from_env()
    t0 = time.time()
    try:
        with psycopg.connect(
                dsn,
                options=f"-c statement_timeout={STATEMENT_TIMEOUT_MS} "
                        f"-c idle_in_transaction_session_timeout={STATEMENT_TIMEOUT_MS} "
                        f"-c search_path=public",
        ) as conn:
            register_vector(conn)
            with conn.cursor() as cur:
                cur.execute(rendered, typed)
                rows = cur.fetchall()
                cols = [d.name for d in cur.description]
    except psycopg.errors.QueryCanceled:
        raise TimeoutError(f"Query timeout (>{STATEMENT_TIMEOUT_MS} ms)")
    except Exception as e:
        raise RuntimeError(f"Database error: {str(e)}")

    elapsed = int((time.time() - t0) * 1000)
    out_rows = [dict(zip(cols, r)) for r in rows]

    return {
        "name": name,
        "columns": cols,
        "rows": out_rows,
        "rowCount": len(out_rows),
        "limit": limit,
        "offset": offset,
        "elapsed_ms": elapsed,
    }


@blp.post("/api/query")
def run_query():
    """
    Run a predefined database query with parameters, paging, and ordering.

    :return: JSON response with query results or error message
    """
    payload = request.get_json(force=True) or {}
    name = payload.get("name")
    params = payload.get("params", {})
    paging = payload.get("paging", {})
    order = payload.get("order", {})

    try:
        result = execute_named_query(name, params=params, paging=paging, order=order)
    except ValueError as e:
        return jsonify({"error": str(e)}), 400
    except TimeoutError as e:
        return jsonify({"error": str(e)}), 408
    except RuntimeError as e:
        return jsonify({"error": str(e)}), 400

    return jsonify(result), 200