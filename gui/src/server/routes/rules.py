"""Read-only access to the default rule set (matching rules + reaction rules), for the GUI's Rules browser tab."""

import threading

import rdkit
from flask import Blueprint, Response, jsonify, request
from rdkit.Chem.Draw import rdMolDraw2D

from retromol.model.rules import ReactionRule, RuleSet

blp_rule_set = Blueprint("rule_set", __name__)

_rule_set: RuleSet | None = None
_reaction_rules_by_id: dict[str, ReactionRule] | None = None

REACTION_SVG_WIDTH = 640
REACTION_SVG_HEIGHT = 220
REACTION_SVG_THEMES = ("light", "dark")

# Rendered reaction scheme SVGs never change for the life of the process (same
# rule set, same fixed canvas size) -- cheap to keep around rather than
# re-running RDKit's layout/drawing on every request for the same (rule, theme).
_reaction_svg_cache: dict[tuple[str, str], str] = {}
_reaction_svg_cache_lock = threading.Lock()


def _get_rule_set() -> RuleSet:
    """
    Lazily load (once per process) and return the default RuleSet.

    Loading reparses the ruleset YAML and rebuilds every rule's RDKit Mol/reaction --
    cheap enough to do once, but wasteful to redo per-request, hence the cache. Kept
    separate from discovery.py's own DiscoveryContext cache, which loads a RuleSet too
    but bundles it with a fingerprint vocabulary and an all-pairs alignment scoring
    matrix -- work this read-only browser endpoint has no reason to force.

    :return: the shared default RuleSet
    """
    global _rule_set
    if _rule_set is None:
        _rule_set = RuleSet.load_default()
    return _rule_set


def _get_reaction_rules_by_id() -> dict[str, ReactionRule]:
    """
    Lazily build (once per process) and return a reaction rule id -> rule lookup.

    :return: the shared id -> ReactionRule map
    """
    global _reaction_rules_by_id
    if _reaction_rules_by_id is None:
        _reaction_rules_by_id = {r.id: r for r in _get_rule_set().reaction_rules}
    return _reaction_rules_by_id


def _render_reaction_svg(rule: ReactionRule, theme: str) -> str:
    """
    Render a reaction rule's scheme (reactant(s) -> product(s)) to an SVG string via
    RDKit, which -- unlike the client-side SMILES-only renderer this replaced --
    understands full SMARTS query syntax (OR-lists like "[C,c]", boolean primitives
    like "[O&D2]"/"[C;!R]") and draws every rule in the set, not just the handful
    whose SMARTS happens to be plain enough for a SMILES grammar.

    :param rule: the reaction rule to render
    :param theme: "light" or "dark" -- selects a monochrome atom/bond color so the
        (background-less) SVG reads correctly against the GUI's current theme
    :return: the rendered SVG markup
    """
    drawer = rdMolDraw2D.MolDraw2DSVG(REACTION_SVG_WIDTH, REACTION_SVG_HEIGHT)
    opts = drawer.drawOptions()
    opts.clearBackground = False
    opts.useBWAtomPalette()
    if theme == "dark":
        opts.setAtomPalette({-1: (1, 1, 1)})
        opts.setSymbolColour((1, 1, 1))
    drawer.DrawReaction(rule.rxn)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def _get_reaction_svg(rule: ReactionRule, theme: str) -> str:
    key = (rule.id, theme)
    cached = _reaction_svg_cache.get(key)
    if cached is not None:
        return cached

    svg = _render_reaction_svg(rule, theme)
    with _reaction_svg_cache_lock:
        _reaction_svg_cache[key] = svg
    return svg


@blp_rule_set.get("/api/ruleSet")
def rule_set() -> tuple[Response, int]:
    """
    The full default rule set (matching rules + reaction rules), for the Rules browser tab.

    Returned in one shot rather than paginated -- the rule set is small (on the order
    of a few hundred entries total) and effectively static for the life of the
    process, same rationale as /api/motifStructures.

    :return: a tuple containing the rule set payload and an HTTP status code
    """
    rules = _get_rule_set()

    matching_rules = [
        {
            "id": r.id,
            "name": r.name,
            "smiles": r.smiles,
            "displaySmiles": r.display_smiles,
            "pseudonyms": r.pseudonyms,
            "terminal": r.terminal,
            "stereochemistry": r.stereochemistry,
            "props": r.props,
        }
        for r in rules.matching_rules
    ]

    reaction_rules = [
        {
            "id": r.id,
            "name": r.name,
            "smarts": r.smarts,
            "allowedInBulk": r.allowed_in_bulk,
            "props": r.props,
        }
        for r in rules.reaction_rules
    ]

    return jsonify({
        "matchStereochemistry": rules.match_stereochemistry,
        "matchingRules": matching_rules,
        "reactionRules": reaction_rules,
    }), 200


@blp_rule_set.get("/api/reactionSchemeSvg/<rule_id>")
def reaction_scheme_svg(rule_id: str) -> tuple[Response, int]:
    """
    A reaction rule's scheme (reactant(s) -> product(s)), rendered server-side by
    RDKit as an SVG -- see _render_reaction_svg for why this replaced the earlier
    client-side (SMILES-only) depiction.

    :param rule_id: the reaction rule's id (ReactionRule.id, a sha256 of its SMARTS)
    :return: a tuple containing {"svg": ...} and an HTTP status code
    """
    theme = request.args.get("theme", "light")
    if theme not in REACTION_SVG_THEMES:
        theme = "light"

    rule = _get_reaction_rules_by_id().get(rule_id)
    if rule is None:
        return jsonify({"error": "Unknown reaction rule id"}), 404

    svg = _get_reaction_svg(rule, theme)
    return jsonify({"svg": svg, "rdkitVersion": rdkit.__version__}), 200
