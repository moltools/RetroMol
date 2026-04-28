"""Visualization utilities for ReactionGraph."""

import re

from retromol.model.reaction_graph import ReactionGraph
from retromol.chem.mol import mol_to_smiles


def visualize_reaction_graph(g: ReactionGraph, html_path: str, root_enc: str | None = None) -> None:
    """
    Visualize ReactionGraph.

    :param g: ReactionGraph to visualize.
    :param html_path: Path to save the HTML visualization.
    :param root_enc: Optional root molecule encoding to highlight.
    .. note:: Requires pyvis package.
    """

    try:
        from pyvis.network import Network
    except ImportError as e:
        raise ImportError("Requires pyvis. Install with: pip install pyvis!") from e

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

        smiles = mol_to_smiles(node.mol, include_tags=False)

        net.add_node(
            mol_vid(enc),
            label=label,
            shape="ellipse",
            color=color,
            smiles=smiles,
        )

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
          "interaction": {
            "navigationButtons": true
          },
          "physics": {"stabilization": true}
        }
        """
    )

    net.write_html(html_path, notebook=False)

    _copy_popup_js = r"""
    // Custom SMILES copy popup
    const smilesPopup = document.createElement("div");
    smilesPopup.id = "smiles-copy-popup";
    smilesPopup.style.position = "fixed";
    smilesPopup.style.display = "none";
    smilesPopup.style.zIndex = "999999";
    smilesPopup.style.background = "white";
    smilesPopup.style.border = "1px solid #ccc";
    smilesPopup.style.borderRadius = "6px";
    smilesPopup.style.padding = "8px";
    smilesPopup.style.boxShadow = "0 2px 8px rgba(0,0,0,0.25)";
    smilesPopup.style.fontFamily = "sans-serif";

    const smilesButton = document.createElement("button");
    smilesButton.textContent = "Copy SMILES";
    smilesButton.style.cursor = "pointer";

    smilesPopup.appendChild(smilesButton);
    document.body.appendChild(smilesPopup);

    let currentSmiles = null;
    let smilesHideTimer = null;

    function showSmilesPopup(x, y, smiles) {
      currentSmiles = smiles;
      smilesButton.textContent = "Copy SMILES";
      smilesPopup.style.left = `${x + 20}px`;
      smilesPopup.style.top = `${y + 20}px`;
      smilesPopup.style.display = "block";
    }

    function hideSmilesPopupSoon() {
      clearTimeout(smilesHideTimer);
      smilesHideTimer = setTimeout(() => {
        smilesPopup.style.display = "none";
        currentSmiles = null;
      }, 700);
    }

    smilesPopup.addEventListener("mouseenter", function () {
      clearTimeout(smilesHideTimer);
    });

    smilesPopup.addEventListener("mouseleave", hideSmilesPopupSoon);

    smilesButton.addEventListener("click", function () {
      if (!currentSmiles) return;

      navigator.clipboard.writeText(currentSmiles).then(
        function () {
          smilesButton.textContent = "Copied!";
        },
        function () {
          // fallback for file:// pages or older browsers
          const ta = document.createElement("textarea");
          ta.value = currentSmiles;
          document.body.appendChild(ta);
          ta.select();
          document.execCommand("copy");
          document.body.removeChild(ta);
          smilesButton.textContent = "Copied!";
        }
      );
    });

    container.addEventListener("mousemove", function (ev) {
      const nodeId = network.getNodeAt({
        x: ev.offsetX,
        y: ev.offsetY,
      });

      if (!nodeId) {
        return;
      }

      const node = nodes.get(nodeId);

      if (!node || !node.smiles) {
        return;
      }

      clearTimeout(smilesHideTimer);
      showSmilesPopup(ev.clientX, ev.clientY, node.smiles);
    });

    container.addEventListener("mouseleave", hideSmilesPopupSoon);
    """

    with open(html_path, "r", encoding="utf-8") as f:
        html = f.read()

    pattern = r"(network\s*=\s*new\s+vis\.Network\s*\(\s*container\s*,\s*data\s*,\s*options\s*\)\s*;)"

    if not re.search(pattern, html):
        raise RuntimeError(
            "Could not find PyVis network creation line. "
            "Open the HTML and search for 'new vis.Network' to inspect the exact line."
        )

    html = re.sub(
        pattern,
        r"\1\n" + _copy_popup_js,
        html,
        count=1,
    )

    with open(html_path, "w", encoding="utf-8") as f:
        f.write(html)
