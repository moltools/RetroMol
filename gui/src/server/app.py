"""Module for configuring the Flask app."""

import logging
import os
import time

from flask import Flask, jsonify
from flask_cors import CORS

from routes.session import (
    blp_create_session,
    blp_delete_session,
    blp_get_session,
    blp_save_session,
    blp_delete_item,
)
from routes.session_store import get_or_init_app_start_epoch
from routes.database import check_database_ready, duckdb_path_from_env
from routes.jobs import blp_search_compound, blp_submit_compound, blp_reconstruct_compound, blp_submit_gene_cluster
from routes.events import blp_events, blp_sse_ticket


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
# In development: simple basicConfig
# In production (under gunicorn): reuse gunicorn's error logger handlers
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

app.logger.info("Flask logger configured")


# Set environment and debug mode
app.config["ENV"] = os.getenv("FLASK_ENV", "production")  # defaults to "production"
app.config["DEBUG"] = app.config["ENV"] == "development"
print("starting app in environment:", app.config["ENV"])
print("debug mode is:", app.debug)


# Log the environment
if app.config["ENV"] == "production":
    print("production environment detected")
elif app.config["ENV"] == "development":
    print("development environment detected")
else:
    print(f"unknown environment: {app.config['ENV']}")


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
    try:
        check_database_ready()
        return jsonify({
            "status": "ready",
            "database": str(duckdb_path_from_env()),
        }), 200
    except Exception as e:
        return jsonify({
            "status": "not ready",
            "error": str(e),
        }), 503


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

app.register_blueprint(blp_events)
app.register_blueprint(blp_sse_ticket)
