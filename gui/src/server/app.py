"""Module for configuring the Flask app."""

import logging
import os
import time

from flask import Flask, jsonify, request
from flask_cors import CORS
from pythonjsonlogger import jsonlogger
from prometheus_flask_exporter import PrometheusMetrics

from routes.session import (
    blp_create_session,
    blp_delete_session,
    blp_get_session,
    blp_save_session,
    blp_delete_item,
)
from routes.session_store import get_or_init_app_start_epoch, redis_client
from routes.database import check_database_ready, duckdb_path_from_env
from routes.stats import blp_database_stats
from routes.jobs import (
    blp_search_compound,
    blp_submit_compound,
    blp_reconstruct_compound,
    blp_submit_gene_cluster,
    blp_get_cluster_readout,
    blp_reconstruct_gene_cluster,
)
from routes.events import blp_events, blp_sse_ticket
from routes.discovery import (
    blp_discovery_monomer_names,
    blp_discovery_query,
    blp_discovery_msa,
    blp_discovery_compare_tanimoto,
    blp_submit_discovery_query,
    blp_get_discovery_query_result,
    blp_motif_structures,
    get_discovery_context,
)
from routes.rate_limit import limiter, RATE_LIMIT_REJECTIONS
from routes.rules import blp_rule_set, blp_generate_backbone
from routes.enrichment import blp_entry_search, blp_enrichment_analysis
from routes.browse import blp_browse_entries, blp_annotation_terms, blp_export_entries


# Initialize the Flask app
app = Flask(__name__)

# Enable CORS in development environment
if os.getenv("FLASK_ENV") == "development":
    origins = ["http://localhost:3000"]
    CORS(
        app,
        resources={r"/api/*": {"origins": origins}},
        supports_credentials=False,  # set True only if you actually use cookies
    )

# Logging setup
# In development: simple basicConfig, human-readable.
# In production (under gunicorn): reuse gunicorn's error logger handlers, but with a
# JSON formatter -- turns "is the concurrency cap actually being hit" from a
# question into something `docker logs | jq` can answer directly.
if os.getenv("FLASK_ENV") == "development":
    logging.basicConfig(level=logging.DEBUG)
    app.logger.setLevel(logging.DEBUG)
else:
    gunicorn_logger = logging.getLogger("gunicorn.error")
    if gunicorn_logger.handlers:
        app.logger.handlers = gunicorn_logger.handlers
        app.logger.setLevel(gunicorn_logger.level)
    else:
        # Fallback if not under gunicorn
        logging.basicConfig(level=logging.INFO)
        app.logger.setLevel(logging.INFO)

    json_formatter = jsonlogger.JsonFormatter(
        "%(asctime)s %(name)s %(levelname)s %(message)s", rename_fields={"asctime": "timestamp", "levelname": "level"}
    )
    for handler in app.logger.handlers:
        handler.setFormatter(json_formatter)

    # Without this, records also propagate to the root logger (which gunicorn
    # configures separately, with its own plain formatter) -- every app.logger call
    # would print twice: once JSON, once plain. app.logger.handlers already covers
    # everywhere these records need to go.
    app.logger.propagate = False

app.logger.info("Flask logger configured")


# Set environment and debug mode
app.config["ENV"] = os.getenv("FLASK_ENV", "production")  # defaults to "production"
app.config["DEBUG"] = app.config["ENV"] == "development"
app.config["START_EPOCH"] = get_or_init_app_start_epoch()
print("starting app in environment:", app.config["ENV"])
print("debug mode is:", app.debug)


# Log the environment
if app.config["ENV"] == "production":
    print("production environment detected")
elif app.config["ENV"] == "development":
    print("development environment detected")
else:
    print(f"unknown environment: {app.config['ENV']}")


# Rate limiting (see routes/rate_limit.py) and Prometheus metrics at /metrics --
# request count/latency/in-progress by endpoint+status are auto-instrumented;
# routes/queue.py and routes/rate_limit.py add a couple of custom counters on top
# for job outcomes and rate-limit rejections specifically.
limiter.init_app(app)
metrics = PrometheusMetrics(app)


@app.after_request
def _log_request(response):
    """
    Log one structured JSON line per request: method, path, status, duration.

    :param response: the outgoing Flask response
    :return: the same response, unmodified
    """
    duration_ms = None
    if hasattr(request, "_start_time"):
        duration_ms = int((time.time() - request._start_time) * 1000)

    app.logger.info(
        "request",
        extra={
            "method": request.method,
            "path": request.path,
            "status": response.status_code,
            "duration_ms": duration_ms,
        },
    )
    return response


@app.before_request
def _mark_request_start():
    """Stamp the request with a start time for _log_request's duration calculation."""
    request._start_time = time.time()


@app.errorhandler(429)
def rate_limited(e) -> tuple[dict, int]:
    """
    Handle rate-limit rejections (Flask-Limiter raises a 429 on breach).

    :param e: the raised exception (its description holds Flask-Limiter's own message)
    :return: a dictionary describing the rejection and an HTTP status code
    """
    RATE_LIMIT_REJECTIONS.labels(path=request.path).inc()
    return jsonify({"error": f"Rate limit exceeded: {e.description}"}), 429


# Register api endpoints
@app.errorhandler(404)
def not_found(_) -> str:
    """
    Handle 404 errors by returning the main index page.

    :param _: the error, not used
    :return: the index HTML page
    """
    return app.send_static_file("index.html")


@app.route("/")
def index() -> str:
    """
    Serve the main index page.

    :return: the index HTML page
    """
    return app.send_static_file("index.html")


@app.route("/api/startup")
def startup() -> tuple[dict[str, int], int]:
    """
    Get the startup time of the server.

    :return: a dictionary with startup, current time, uptime and HTTP status code
    """
    startup_epoch = get_or_init_app_start_epoch()
    return jsonify({
        "startup": startup_epoch,
        "current": int(time.time()),
        "uptime": int(time.time()) - startup_epoch,
    }), 200


@app.route("/api/health", methods=["GET"])
def health() -> tuple[dict[str, str], int]:
    """
    Health check endpoint.

    :return: a dictionary indicating the server is healthy and HTTP status code
    """
    return jsonify({
        "status": "ok",
        "time": int(time.time()),
        "uptime": int(time.time()) - app.config["START_EPOCH"],
    }), 200


@app.route("/api/ready", methods=["GET"])
def ready() -> tuple[dict[str, str], int]:
    """
    Readiness check: confirms both DuckDB and Redis are actually reachable, not just
    that the process is up -- what container healthchecks are wired to (see
    gui/docker/backend.Dockerfile).

    :return: a dictionary indicating readiness (or which dependency failed) and an HTTP status code
    """
    try:
        check_database_ready()
    except Exception as e:
        return jsonify({"status": "not ready", "error": f"duckdb: {e}"}), 503

    try:
        redis_client.ping()
    except Exception as e:
        return jsonify({"status": "not ready", "error": f"redis: {e}"}), 503

    return jsonify({
        "status": "ready",
        "database": str(duckdb_path_from_env()),
    }), 200


# Register blueprints
app.register_blueprint(blp_create_session)
app.register_blueprint(blp_delete_session)
app.register_blueprint(blp_get_session)
app.register_blueprint(blp_save_session)
app.register_blueprint(blp_delete_item)
app.register_blueprint(blp_search_compound)
app.register_blueprint(blp_submit_compound)
app.register_blueprint(blp_reconstruct_compound)
app.register_blueprint(blp_submit_gene_cluster)
app.register_blueprint(blp_get_cluster_readout)
app.register_blueprint(blp_reconstruct_gene_cluster)
app.register_blueprint(blp_database_stats)

app.register_blueprint(blp_events)
app.register_blueprint(blp_sse_ticket)
app.register_blueprint(blp_discovery_monomer_names)
app.register_blueprint(blp_motif_structures)
app.register_blueprint(blp_discovery_query)
app.register_blueprint(blp_discovery_msa)
app.register_blueprint(blp_discovery_compare_tanimoto)
app.register_blueprint(blp_submit_discovery_query)
app.register_blueprint(blp_get_discovery_query_result)
app.register_blueprint(blp_rule_set)
app.register_blueprint(blp_generate_backbone)
app.register_blueprint(blp_entry_search)
app.register_blueprint(blp_enrichment_analysis)
app.register_blueprint(blp_browse_entries)
app.register_blueprint(blp_annotation_terms)
app.register_blueprint(blp_export_entries)

# The two rate-limit tiers on top of the app-wide default (see routes/rate_limit.py)
# are applied as @limiter.limit(...) decorators directly on each route function, in
# jobs.py/discovery.py -- Flask-Limiter's documented pattern for a limiter
# instance created in a separate module from the app (decorate at definition time,
# call limiter.init_app(app) once the app exists, as done above).

# Warm the discovery context cache eagerly so the first user request isn't slow --
# building it reparses the ruleset and computes an all-pairs Tanimoto scoring matrix
# over every named monomer. Non-fatal on failure: it'll just build lazily on first use.
try:
    get_discovery_context()
    app.logger.info("Discovery context warmed successfully")
except Exception:
    app.logger.exception("Failed to warm discovery context at startup; will build lazily on first request")
