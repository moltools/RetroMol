"""Module for RetroMol results readout."""

from typing import Any

from networkx import (
    Graph,
    all_pairs_shortest_path_length,
    connected_components,
    is_connected,
    shortest_path,
)

from retromol.api import timeout_decorator
from retromol.chem import Mol, smiles_to_mol
from retromol.config import TIMEOUT_LINEAR_READOUT, TIMEOUT_OPTIMAL_MAPPINGS
from retromol.graph import merge_nodes, mol_to_graph
from retromol.io import Result


def optimal_mappings(result: Result) -> list[dict[str, Any]]:
    """
    Return one mapping per nesting level.
    At each level, nodes in the new level take precedence.
    Previous level nodes are added only if they do not overlap
    and are not submappings of the new level nodes.

    :param result: RetroMol Result object
    :return: list of mapping dictionaries, each with keys:
        - "nodes": list of nodes, each with "identity", "smiles", and "tags"
        - "covered_tags": list of all tags covered by the selected nodes
        - "n_nodes": number of nodes in the mapping
        - "n_tags": number of unique tags covered
    """
    identified: set[tuple[str, str, tuple[int, ...]]] = result.get_identified_nodes()

    # Normalize and filter empties
    items: list[dict[str, Any]] = []
    for identity, smiles, tags in identified:
        tag_set = frozenset(int(t) for t in tags)
        if not tag_set:
            continue
        items.append({"identity": identity, "smiles": smiles, "tags": tag_set})

    if not items:
        return []

    # Sort: largest pieces first
    items.sort(key=lambda m: (-len(m["tags"]), str(m["smiles"]), str(m["identity"])))

    # Partition into levels (graph coloring by overlap)
    levels: list[list[dict[str, Any]]] = []
    for item in items:
        placed = False
        for lvl in levels:
            if any(item["tags"] & other["tags"] for other in lvl):
                continue
            lvl.append(item)
            placed = True
            break
        if not placed:
            levels.append([item])

    results: list[dict[str, Any]] = []

    # Process levels with precedence for newer levels
    for lvl_idx, lvl in enumerate(levels):
        used_tags: set[int] = set()
        chosen: list[dict[str, Any]] = []

        # Step 1: take all new-level items
        for m in lvl:
            chosen.append(m)
            used_tags |= m["tags"]

        # Step 2: consider all previous levels
        for prev_idx in range(lvl_idx):
            for m in levels[prev_idx]:
                # skip if overlaps with any new-level item
                if m["tags"] & used_tags:
                    continue
                # skip if fully contained in a new-level item (submapping)
                if any(m["tags"] <= new_m["tags"] for new_m in lvl):
                    continue
                chosen.append(m)
                used_tags |= m["tags"]

        # Build solution dict
        covered: set[int] = set().union(*(n["tags"] for n in chosen))
        nodes = [{"identity": n["identity"], "smiles": n["smiles"], "tags": sorted(n["tags"])} for n in chosen]
        results.append(
            {
                "nodes": nodes,
                "covered_tags": sorted(covered),
                "n_nodes": len(nodes),
                "n_tags": len(covered),
            }
        )

    return results


def mapping_to_graph(tagged_smi: str, mapping: dict[str, Any]) -> "Graph[str | int]":
    """
    Convert a mapping dictionary to a NetworkX graph.

    :param tagged_smi: the SMILES representation of the molecule
    :param mapping: a dictionary representing a mapping with keys:
        - "nodes": list of nodes, each with "identity", "smiles", and "tags"
        - "covered_tags": list of all tags covered by the selected nodes
        - "n_nodes": number of nodes in the mapping
        - "n_tags": number of unique tags covered
    :return: a NetworkX graph representing the mapping
    """
    tagged_mol: Mol = smiles_to_mol(tagged_smi)
    G = mol_to_graph(tagged_mol, use_tags=True)
    G.graph["tagged_smiles"] = tagged_smi

    # Merge nodes in graph
    for node_idx, node in enumerate(mapping["nodes"]):
        node_identity = node["identity"]
        node_smiles = node["smiles"]
        node_tags = node["tags"]
        merge_nodes(
            G,
            merged_node_id=f"node_{node_idx}_{node_identity.replace(' ', '_')}",
            nodes=node_tags,
            props={"identity": node_identity, "smiles": node_smiles, "tags": node_tags},
        )

    # Delete unmerged nodes
    to_remove = [n for n, d in G.nodes(data=True) if "identity" not in d]
    G.remove_nodes_from(to_remove)

    return G


def _dfs_graphs(root: "Graph[int | str]") -> "list[Graph[int | str]]":
    """
    Return all graphs in DFS discovery order via node attr 'graph'.

    :param root: root graph
    :return: list of graphs in DFS order
    """
    out: list[Graph[int | str]] = []

    def dfs(G: "Graph[int | str]"):
        out.append(G)
        # Stable order over nodes for determinism
        for _, d in list(G.nodes(data=True)):
            sub = d.get("graph")
            if isinstance(sub, Graph):
                dfs(sub)

    dfs(root)
    return out


def _bfs_graphs_with_depth(root: "Graph[int | str]") -> "list[tuple[Graph[int | str], int]]":
    """
    Return all graphs with their true nesting depth (root depth=0) using BFS.

    :param root: root graph
    :return: list of tuples (graph, depth)
    """
    out: list[tuple[Graph[int | str], int]] = []
    seen: set[int] = set()
    q: list[tuple[Graph[int | str], int]] = [(root, 0)]
    while q:
        G, depth = q.pop(0)
        if id(G) in seen:
            continue
        seen.add(id(G))
        out.append((G, depth))
        for _, d in G.nodes(data=True):
            sub = d.get("graph")
            if isinstance(sub, Graph):
                q.append((sub, depth + 1))
    return out


def _graphs_with_metadata(root: "Graph[int | str] | None") -> list[dict[str, Any]]:
    """
    Combine DFS order (legacy 'level_index' semantics) with true depth.

    :param root: root graph
    :return: a list of dicts: {"graph": G, "dfs_index": i, "depth": d}
    """
    if root is None:
        return []
    dfs_list = _dfs_graphs(root)
    bfs_list = _bfs_graphs_with_depth(root)  # (G, depth)
    depth_by_id = {id(G): depth for (G, depth) in bfs_list}
    meta: list[dict[str, Any]] = []
    for i, G in enumerate(dfs_list):
        meta.append({"graph": G, "dfs_index": i, "depth": depth_by_id.get(id(G), 0)})
    return meta


def list_levels_summary(result: Result) -> list[dict[str, Any]]:
    """
    Legacy-style summary (DFS order), now also shows true nesting depth per entry.

    :param result: RetroMol Result object.
    :return: ;ist of dicts, each with keys:
        - "dfs_index": DFS discovery index
        - "depth": True nesting depth (root=0)
        - "n_nodes": number of nodes in the graph at this level
        - "n_identified": number of identified nodes (with 'identity' attribute)
    """
    # NOTE: should never happen, but be defensive
    if result.graph is None:
        return []
    summary: list[dict[str, Any]] = []
    for m in _graphs_with_metadata(result.graph):
        G = m["graph"]
        summary.append(
            {
                "dfs_index": m["dfs_index"],
                "depth": m["depth"],
                "n_nodes": G.number_of_nodes(),
                "n_identified": sum(1 for _, d in G.nodes(data=True) if d.get("identity") is not None),
                "wave_names": sorted({d.get("wave_name") for _, d in G.nodes(data=True) if "wave_name" in d}),
            }
        )
    return summary


def find_depth_by_wave_name(result: Result, wave_name: str) -> int | None:
    """
    Return the smallest nesting depth that contains any node with the given wave_name.

    :param result: RetroMol Result object
    :param wave_name: the wave_name to search for
    :return: the smallest nesting depth containing the wave_name, or None if not found
    """
    # NOTE: should never happen, but be defensive
    if result.graph is None:
        return None
    depths: list[int] = []
    for m in _graphs_with_metadata(result.graph):
        G = m["graph"]
        if any(d.get("wave_name") == wave_name for _, d in G.nodes(data=True)):
            depths.append(m["depth"])
    return min(depths) if depths else None


def _monomer_nodes_at_level(G: "Graph[int | str]", require_identified: bool) -> list[Any]:
    """
    Return monomer nodes at the given graph level.

    :param G: the graph at the current level
    :param require_identified: if True, only include nodes with an 'identity' attribute
    :return: list of monomer node identifiers
    """
    nodes: list[Any] = []
    for n, d in G.nodes(data=True):
        if require_identified:
            if d.get("identity") is not None and "tags" in d:
                nodes.append(n)
        else:
            if "tags" in d:
                nodes.append(n)
    return nodes


def _is_path_component(H: "Graph[int | str]", nodes: list[Any]) -> bool:
    """
    Check if the subgraph induced by 'nodes' in H is a path.

    :param H: the host graph
    :param nodes: list of nodes to check
    :return: True if the subgraph is a path, False otherwise
    """
    if not nodes:
        return False
    C = H.subgraph(nodes)
    if not is_connected(C):
        return False
    degs: dict[str | int, int] = dict(C.degree())
    if any(d > 2 for d in degs.values()):
        return False
    if len(C) == 1:
        return True
    ones = sum(1 for d in degs.values() if d == 1)
    return (C.number_of_edges() == len(C) - 1) and (ones == 2)


def _order_nodes_along_path(H: "Graph[int | str]", nodes: list[Any]) -> list[Any]:
    """
    Order nodes along a path component.

    :param H: the host graph
    :param nodes: list of nodes in the path component
    :return: ordered list of nodes along the path
    """
    C = H.subgraph(nodes).copy()
    if len(C) == 1:
        return list(C.nodes())
    degs: dict[str | int, int] = dict(C.degree())
    start = [n for n, d in degs.items() if d == 1][0]
    order = [start]
    prev = None
    cur = start
    while True:
        nbrs = [v for v in C.neighbors(cur) if v != prev]
        if not nbrs:
            break
        nxt = nbrs[0]
        order.append(nxt)
        prev, cur = cur, nxt
    return order


def _longest_path_approx(H: "Graph[int | str]", nodes: list[Any]) -> list[Any]:
    """
    Approximate longest path in the subgraph induced by 'nodes' in H.

    :param H: the host graph
    :param nodes: list of nodes to consider
    :return: list of nodes along the approximate longest path
    """
    C: Graph[int | str] = H.subgraph(nodes).copy()
    if len(C) <= 1:
        return list(C.nodes())
    lengths = dict(all_pairs_shortest_path_length(C))
    max_d = -1
    pair: tuple[Any, Any] | None = None
    for u, dmap in lengths.items():
        for v, d in dmap.items():
            if d > max_d:
                max_d = d
                pair = (u, v)
    assert pair is not None, "At least one pair should exist in non-empty connected graph"
    return shortest_path(C, source=pair[0], target=pair[1])


def _payload_from_order(G_src: "Graph[int | str]", order: list[Any]) -> dict[str, Any]:
    """
    Create payload dictionary from ordered nodes.

    :param G_src: source graph
    :param order: ordered list of nodes
    :return: payload dictionary
    """
    items: list[dict[str, Any]] = []
    for n in order:
        d = G_src.nodes[n]
        items.append(
            {
                "node": n,
                "identity": d.get("identity"),
                "smiles": d.get("smiles"),
                "tags": sorted(d.get("tags", [])),
            }
        )
    return {
        "n_monomers": len(items),
        "ordered_monomers": items,
    }


def _score_payload(pl: dict[str, Any]) -> tuple[int, int, tuple[str, ...]]:
    """
    Score payload for comparison: (n_monomers, n_identified, node_key).

    :param pl: payload dictionary
    :return: scoring tuple
    """
    n = pl["n_monomers"]
    n_ident = sum(1 for it in pl["ordered_monomers"] if it.get("identity") is not None)
    node_key = tuple(str(it["node"]) for it in pl["ordered_monomers"])
    return (n, n_ident, node_key)


def linear_readout(
    result: Result,
    require_identified: bool = True,
    mode: str = "all",  # "all" | "best_per_level" | "global_best"
    nesting_depth: int | None = None,
) -> dict[str, Any]:
    """
    Linear backbone readouts.

    :param result: RetroMol Result object
    :param require_identified: if ``True``, only consider monomer nodes with an
        assigned identity. If ``False``, consider all monomer nodes
    :param nesting_depth: maximum nesting level to analyze. If ``None``, all depths are included
        - when ``nesting_depth`` is ``None``: iterate all graphs in DFS order and return
        a structure identical to the previous version (keys and shapes), but each
        entry now also includes a ``depth`` field for clarity.
        - set ``nesting_depth = k`` to restrict analysis to graphs at that **true**
        nesting level (root = 0, its children = 1, etc.).
    :param mode: determines the aggregation mode of readouts
        Supported values:
        - ``"all"``: return all depth levels and paths.
        - ``"best_per_level"``: return the best backbone per depth level.
        - ``"global_best"``: return the globally best backbone only.

    :returns:
        Depending on the selected ``mode``:

        **mode = "all"**
            Returns:
            ``{"levels": [
                {"dfs_index": int, "depth": int,
                "strict_paths": [payload, ...],
                "fallback": payload_or_None},
                ...
            ]}``

        **mode = "best_per_level"**
            Returns:
            ``{"levels": [
                {"dfs_index": int, "depth": int,
                "strict_path": bool,
                "backbone": payload,
                "notes": str},
                ...
            ]}``

        **mode = "global_best"**
            Returns:
            ``{"dfs_index": int, "depth": int,
            "strict_path": bool,
            "backbone": payload,
            "notes": str}``
    :rtype:
        Dict[str, Any]
    """
    metas = _graphs_with_metadata(result.graph)
    if nesting_depth is not None:
        metas = [m for m in metas if m["depth"] == nesting_depth]
        if not metas:
            msg = f"No graphs at nesting_depth={nesting_depth}."
            if mode == "global_best":
                return {
                    "dfs_index": -1,
                    "depth": nesting_depth,
                    "strict_path": False,
                    "backbone": {"n_monomers": 0, "ordered_monomers": []},
                    "notes": msg,
                }
            else:
                return {"levels": [], "notes": msg}

    # Per-graph analysis
    entries: list[dict[str, Any]] = []
    for m in metas:
        G = m["graph"]
        dfs_idx = m["dfs_index"]
        depth = m["depth"]

        monomer_nodes = _monomer_nodes_at_level(G, require_identified)
        if not monomer_nodes:
            entries.append(
                {
                    "dfs_index": dfs_idx,
                    "depth": depth,
                    "strict_paths": [],
                    "fallback": None,
                }
            )
            continue

        MG = G.subgraph(monomer_nodes).copy()
        comps = list(connected_components(MG))

        strict_payloads: list[dict[str, Any]] = []
        for comp in comps:
            nodes = list(comp)
            if _is_path_component(MG, nodes):
                order = _order_nodes_along_path(MG, nodes)
                strict_payloads.append(_payload_from_order(G, order))

        fallback_payload = None
        if not strict_payloads and comps:
            largest = max(comps, key=len)
            approx_order = _longest_path_approx(MG, list(largest))
            fallback_payload = _payload_from_order(G, approx_order)

        entries.append(
            {
                "dfs_index": dfs_idx,
                "depth": depth,
                "strict_paths": strict_payloads,
                "fallback": fallback_payload,
            }
        )

    # Assemble per mode
    if mode == "all":
        # Preserve old shape (list under "levels"), but each entry includes depth now.
        # Keep DFS ordering for stability.
        entries.sort(key=lambda e: e["dfs_index"])
        return {"levels": entries}

    if mode == "best_per_level":
        best_levels: list[dict[str, Any]] = []
        for e in sorted(entries, key=lambda x: x["dfs_index"]):
            candidates: list[tuple[bool, dict[str, Any]]] = [(True, pl) for pl in e["strict_paths"]]
            if not candidates and e["fallback"] is not None:
                candidates.append((False, e["fallback"]))
            if not candidates:
                continue
            best_idx = max(range(len(candidates)), key=lambda i: _score_payload(candidates[i][1]))
            is_strict, payload = candidates[best_idx]
            best_levels.append(
                {
                    "dfs_index": e["dfs_index"],
                    "depth": e["depth"],
                    "strict_path": is_strict,
                    "backbone": payload,
                    "notes": "Strict path" if is_strict else "Fallback to longest-path approximation",
                }
            )
        return {"levels": best_levels}

    if mode == "global_best":
        best: tuple[bool, dict[str, Any], int, int] | None = None  # (is_strict, payload, dfs_index, depth)
        for e in entries:
            for pl in e["strict_paths"]:
                if best is None or _score_payload(pl) > _score_payload(best[1]):
                    best = (True, pl, e["dfs_index"], e["depth"])
            if not e["strict_paths"] and e["fallback"] is not None:
                pl = e["fallback"]
                if best is None or _score_payload(pl) > _score_payload(best[1]):
                    best = (False, pl, e["dfs_index"], e["depth"])
        if best is None:
            return {
                "dfs_index": -1,
                "depth": nesting_depth if nesting_depth is not None else -1,
                "strict_path": False,
                "backbone": {"n_monomers": 0, "ordered_monomers": []},
                "notes": "No monomer backbones found.",
            }
        is_strict, pl, dfs_idx, depth = best
        return {
            "dfs_index": dfs_idx,
            "depth": depth,
            "strict_path": is_strict,
            "backbone": pl,
            "notes": "Strict path" if is_strict else "Fallback to longest-path approximation",
        }

    raise ValueError(f"Unknown mode: {mode!r}")


# Decorate the optimal_mappings function with a timeout
optimal_mappings_with_timeout = timeout_decorator(seconds=TIMEOUT_OPTIMAL_MAPPINGS)(optimal_mappings)


# Decorate the linear_readout function with a timeout
linear_readout_with_timeout = timeout_decorator(seconds=TIMEOUT_LINEAR_READOUT)(linear_readout)
