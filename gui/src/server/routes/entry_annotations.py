"""Entry annotations: every phylogeny/biosynthetic_class/chemical_class/bioactivity
term linked to one entry, each with a link to that term's source website where one
exists -- used by the Discovery tab's expanded result view ("what do we know about
this compound/BGC").
"""

from flask import Blueprint, Response, jsonify, request

from routes.database import open_retromol_db

blp_entry_annotations = Blueprint("entry_annotations", __name__)

NCBI_TAXONOMY_URL = "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id={id}"
CHEBI_ENTITY_URL = "https://www.ebi.ac.uk/chebi/searchId.do?chebiId={id}"


def _annotation_url(category: str, rank: str | None, external_id: str | None) -> str | None:
    """Build a "view on <source>" link for one term, if its category/rank has a known
    external database and it resolved an id there. biosynthetic_class and chemical_class
    (NPClassifier) have no canonical per-label public page, so those always return None."""
    if not external_id:
        return None
    if category == "phylogeny":
        return NCBI_TAXONOMY_URL.format(id=external_id)
    if category == "bioactivity":
        if rank in ("chebi_biological_role", "chebi_chemical_role"):
            return CHEBI_ENTITY_URL.format(id=external_id)
    return None


@blp_entry_annotations.get("/api/entryAnnotations")
def entry_annotations() -> tuple[Response, int]:
    """
    Every annotation term linked to one entry, for the Discovery tab's expanded result view.

    :return: a tuple containing the annotation list and an HTTP status code
    """
    entry_id = (request.args.get("entryId") or "").strip()
    if not entry_id:
        return jsonify({"error": "entryId is required"}), 400

    try:
        with open_retromol_db() as db:
            terms = db.entry_annotation_terms(entry_id)
    except Exception as e:
        return jsonify({"error": str(e)}), 503

    return (
        jsonify(
            {
                "results": [
                    {
                        "id": t.id,
                        "category": t.category,
                        "rank": t.rank,
                        "label": t.label,
                        "externalId": t.external_id,
                        "url": _annotation_url(t.category, t.rank, t.external_id),
                    }
                    for t in terms
                ]
            }
        ),
        200,
    )
