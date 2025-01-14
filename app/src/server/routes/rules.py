# -*- coding: utf-8 -*-

"""API endpoints for submitting jobs."""

from flask import Blueprint, Response

from retromol.data.rules import LINEARIZATION_RULES, SEQUENCING_RULES

from .common import ResponseData, Status


# Blueprint for submitting a compound job
blueprint_get_reaction_rules = Blueprint("get_reaction_rules", __name__)


@blueprint_get_reaction_rules.route("/api/get_reaction_rules", methods=["POST"])
def get_reaction_rules() -> Response:
    """Return the reaction rules."""
    payload={
        "preprocessing_rules": LINEARIZATION_RULES,
        "sequencing_rules": SEQUENCING_RULES
    }
    return ResponseData(status=Status.SUCCESS, payload=payload, message="Reaction rules successfully retrieved.").to_dict()