# -*- coding: utf-8 -*-

"""API endpoints for submitting jobs."""

from flask import Blueprint, Response, request, redirect

from retromol.api import run_retromol

from .common import ResponseData, Status


# Blueprint for submitting a compound job
blueprint_submit_compound = Blueprint("submit_compound", __name__)


@blueprint_submit_compound.route("/api/submit_compound", methods=["POST"])
def submit_compound() -> Response:
    """Submit a compound job."""
    # Get the data
    data = request.get_json().get("data", None)

    if data is None:
        return ResponseData(status=Status.FAILURE, message="No data provided.").to_dict()

    # Check if the data is valid
    if "smiles" not in data:
        return ResponseData(status=Status.FAILURE, message="No SMILES provided.").to_dict()

    else:
        smiles = data["smiles"]
        payload = run_retromol("compound", smiles)
        return ResponseData(status=Status.SUCCESS, payload=payload, message="Compound parsed.").to_dict()