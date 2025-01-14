# -*- coding: utf-8 -*-

"""API for RetroMol."""

from flask import Response
from routes.app import app
from routes.submit import blueprint_submit_compound
from routes.rules import blueprint_get_reaction_rules


# Register routes
app.register_blueprint(blueprint_submit_compound)
app.register_blueprint(blueprint_get_reaction_rules)


@app.errorhandler(404)
def not_found(error) -> Response:
    return app.send_static_file("index.html")


@app.route("/")
def index() -> Response:
    return app.send_static_file("index.html")


# api endpoint for fetching version
@app.route("/api/version")
def version() -> Response:
    return {"version": "1.0.0-dev"}
