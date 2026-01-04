"""Visualization utilities for ReactionGraph."""

from retromol.model.reaction_graph import ReactionGraph
from retromol.chem.mol import mol_to_smiles


def visualize_reaction_graph(g: ReactionGraph, html_path: str, root_enc: str | None = None) -> None:
    """
    Visualize ReactionGraph.

    :param g: ReactionGraph to visualize
    :param html_path: path to save the HTML visualization
    :param root_enc: optional root molecule encoding to highlight
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
    def mol_vid(enc: str) -> str:
        """
        Generate a unique molecule node ID.
        
        :param enc: molecule encoding
        :return: str: unique molecule node ID
        """
        return f"m:{enc}"

    def rxn_vid(i: int) -> str:
        """
        Generate a unique reaction node ID.
        
        :param i: reaction index
        :return: str: unique reaction node ID
        """
        return f"r:{i}"

    # Add molecule nodes
    for enc, node in g.nodes.items():

        color = "lightgreen" if root_enc is not None and enc == root_enc else "lightblue"

        identity = None
        if enc not in identified and enc == root_enc:
            identity = "root"
        elif enc in identified and identified[enc]:
            identity = identified[enc].name

        label = str(identity) if identity else "mol"
        net.add_node(mol_vid(enc), label=label, title=label, shape="ellipse", color=color, smiles=mol_to_smiles(node.mol, include_tags=False))

    # Add reaction nodes, and edges between molecules and reactions
    for i, e in enumerate(g.edges):
        title = ", ".join(e.step.names) if getattr(e.step, "names", None) else ""

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
