"""Visualization utilities for ReactionGraph."""

from retromol.model.reaction_graph import ReactionGraph


def visualize_reaction_graph(g: ReactionGraph, html_path: str) -> None:
    """
    Visualize ReactionGraph.

    :param g: ReactionGraph to visualize
    :param html_path: path to save the HTML visualization
    .. note:: requires pyvis package
    """
    try:
        from pyvis.network import Network
    except ImportError as e:
        raise ImportError("Requires pyvis. Install with: pip install pyvis") from e

    # Build identified map from your graph (as in your code)
    identified = {}
    for enc, node in getattr(g, "identified_nodes", {}).items():
        identified[enc] = node.identity

    net = Network(height="800px", width="100%", directed=True, notebook=False)
    net.toggle_physics(True)

    # Use prefixed IDs to avoid collisions with reaction node IDs
    def mol_vid(enc: int) -> str:
        return f"m:{enc}"

    def rxn_vid(i: int) -> str:
        return f"r:{i}"

    # Add molecule nodes
    for enc, node in g.nodes.items():

        identity = None
        if enc in identified and identified[enc]:
            identity = identified[enc].get("identity", None)

        label = str(identity) if identity else "mol"
        net.add_node(mol_vid(enc), label=label, title=label, shape="ellipse")

    # Add reaction nodes, and edges between molecules and reactions
    for i, e in enumerate(g.edges):
        title = ", ".join(e.step.rids) if getattr(e.step, "rids", None) else ""

        net.add_node(rxn_vid(i), label="rxn", title=title, shape="box")

        # src mol -> reaction
        if e.src in g.nodes:
            net.add_edge(mol_vid(e.src), rxn_vid(i), title="reactant", arrows="to")

        # reaction -> dst mol(s)
        for dst in e.dsts:
            if dst not in g.nodes:
                continue
            net.add_edge(rxn_vid(i), mol_vid(dst), title="product", arrows="to")

    # Options
    net.set_options(
        """
        var options = {
          "edges": {"smooth": false},
          "interaction": {"hover": true, "tooltipDelay": 80},
          "physics": {"stabilization": true}
        }
        """
    )

    net.write_html(html_path, notebook=False)
