"""Module contains utilities for defining and working with assembly graphs."""

from dataclasses import dataclass, asdict
from typing import Any, Iterable, Iterable, Iterator, Generator

from rdkit.Chem.rdchem import Mol
import matplotlib.pyplot as plt
import networkx as nx

from retromol.model.reaction_graph import MolNode
from retromol.model.identity import MolIdentity
from retromol.chem.tagging import get_tags_mol


@dataclass(frozen=True)
class RootBondLink:
    """
    One root bond connects two monomer tag-sets.
    
    :var a1_idx: index of the first atom in the root bond
    :var a2_idx: index of the second atom in the root bond
    :var a1_tag: tag of the first atom in the root bond
    :var a2_tag: tag of the second atom in the root bond
    :var a1_symbol: element symbol of the first atom
    :var a2_symbol: element symbol of the second atom
    :var bond_type: stringified version of RDKit BondType
    :var bond_order: bond order if available
    """

    a1_idx: int
    a2_idx: int
    a1_tag: int
    a2_tag: int

    a1_symbol: str
    a2_symbol: str

    bond_type: str  # stringified version of RDKit BondType
    bond_order: float | int | None  # include if available

    def to_dict(self) -> dict[str, Any]:
        """
        Convert the RootBondLink to a dictionary.

        :return: dictionary representation of the RootBondLink
        """
        return asdict(self)
    
    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "RootBondLink":
        """
        Create a RootBondLink from a dictionary.

        :param data: dictionary representation of the RootBondLink
        :return: RootBondLink instance
        """
        return cls(**data)


def build_assembly_graph(
    root_mol: Mol,
    monomers: Iterable[MolNode],
    allow_overlaps: bool = False,
    include_unassigned: bool =False,
) -> nx.Graph:
    """
    Build an assembly graph from the given root molecule and monomers.

    :param root_mol: RDKit Mol representing the root molecule
    :param monomers: iterable of MolNode representing the monomers
    :param allow_overlaps: whether to allow overlapping monomers (default: False)
    :param include_unassigned: whether to include unassigned regions as a node (default: False)
    :return: NetworkX graph representing the assembly graph
    """
    monomers = list(monomers)

    tag_to_monomer: dict[int, str] = {}
    monomer_to_tags: dict[str, set[int]] = {}

    # Map root tags -> monomers
    for m in monomers:
        tags = get_tags_mol(m.mol)
        monomer_to_tags[m.enc] = tags

        for t in tags:
            if t in tag_to_monomer and not allow_overlaps:
                raise ValueError(f"root tag {t} appears in multiple monomersL {tag_to_monomer[t]} and {m.enc}")
            tag_to_monomer[t] = m.enc

    # Initialize empty graph
    g = nx.Graph()

    # Add monomer nodes
    for m in monomers:
        tags = monomer_to_tags[m.enc]
        identity = m.identity
        g.add_node(m.enc, molnode=m, tags=tags, identity=identity)

    UNASSIGNED = "unassigned"
    if include_unassigned:
        g.add_node(UNASSIGNED, molnode=None, tags=set(), identity=None)

    # Scan root bonds
    for b in root_mol.GetBonds():
        a1 = b.GetBeginAtom()
        a2 = b.GetEndAtom()

        t1 = int(a1.GetIsotope())
        t2 = int(a2.GetIsotope())

        if t1 == 0 or t2 == 0:
            continue  # skip non-tagged atoms

        m1 = tag_to_monomer.get(t1)
        m2 = tag_to_monomer.get(t2)

        if m1 is None or m2 is None:
            if not include_unassigned:
                continue  # skip bonds to unassigned regions
            m1 = m1 if m1 is not None else UNASSIGNED
            m2 = m2 if m2 is not None else UNASSIGNED

        if m1 == m2:
            continue  # skip intra-monomer bonds

        link = RootBondLink(
            a1_idx=a1.GetIdx(),
            a2_idx=a2.GetIdx(),
            a1_tag=t1,
            a2_tag=t2,
            a1_symbol=a1.GetSymbol(),
            a2_symbol=a2.GetSymbol(),
            bond_type=str(b.GetBondType()),
            bond_order=float(b.GetBondTypeAsDouble()) if hasattr(b, "GetBondTypeAsDouble") else None,
        )

        if g.has_edge(m1, m2):
            g[m1][m2]["bonds"].append(link)
            g[m1][m2]["n_bonds"] += 1
        else:
            g.add_edge(m1, m2, bonds=[link], n_bonds=1)

    return g


@dataclass(frozen=True, slots=True)
class AssemblyGraph:
    """
    Assembly graph representing monomer connectivity in a molecule.

    :var g: NetworkX graph representing the assembly graph
    :var unassigned: name of the unassigned node
    :var validate: validate graph structure upon initialization
    """

    g: nx.Graph
    unassigned: str = "unassigned"
    validate_upon_initialization: bool = False

    def __post_init__(self) -> None:
        """
        Post-initialization to validate the graph if requested.
        """
        if self.validate_upon_initialization:
            self.validate()

    def __str__(self) -> str:
        """
        String representation of the AssemblyGraph.
        
        :return: string representation
        """
        return f"AssemblyGraph(num_nodes={self.g.number_of_nodes()}, num_edges={self.g.number_of_edges()})"
    
    def monomer_ids(self) -> list[str]:
        """
        Get the list of monomer IDs in the assembly graph.

        :return: list of monomer IDs
        """
        return [n for n in self.g.nodes if n != self.unassigned]
    
    def monomer_nodes(self) -> list[MolNode]:
        """
        Get the list of monomer MolNodes in the assembly graph.

        :return: list of MolNode instances
        """
        out: list[MolNode] = []
        for n in self.monomer_ids():
            mn = self.g.nodes[n]["molnode"]
            if mn is None:
                raise ValueError(f"AssemblyGraph node {n!r} has None molnode")
            out.append(mn)

        return out
    
    def edges_with_bonds(self) -> Iterator[tuple[str, str, list[RootBondLink]]]:
        """
        Iterate over edges with their associated root bonds.

        :return: iterator of tuples (node1, node2, list of RootBondLink)
        """
        for u, v, data in self.g.edges(data=True):
            yield u, v, data["bonds"]

    def drop_unassigned(self) -> "AssemblyGraph":
        """
        Drop the unassigned node from the assembly graph.
        
        :return: AssemblyGraph without the unassigned node
        """
        h = self.g.copy()

        if h.has_node(self.unassigned):
            h.remove_node(self.unassigned)

        return AssemblyGraph(g=h, unassigned=self.unassigned, validate_upon_initialization=True)
    
    def drop_singletons(self) -> "AssemblyGraph":
        """
        Drop singleton nodes (nodes with degree 0) from the assembly graph.
        
        :return: AssemblyGraph with singleton nodes removed
        """
        h = self.g.copy()

        singletons = [n for n, d in h.degree() if d == 0]
        h.remove_nodes_from(singletons)

        return AssemblyGraph(g=h, unassigned=self.unassigned, validate_upon_initialization=True)
    
    def filtered_by_root_bond_elements(
        self,
        allow_pairs: set[frozenset[str]] | None = None,
        drop_isolated: bool = True,
    ) -> "AssemblyGraph":
        """
        Filter the assembly graph by allowed root bond element pairs.
        
        :param allow_pairs: set of allowed element symbol pairs (as frozensets)
        :param drop_isolated: whether to drop isolated nodes after filtering (default: True)
        :return: filtered AssemblyGraph
        """
        h = self.g.copy()

        if allow_pairs is None:
            return AssemblyGraph(g=h, unassigned=self.unassigned)
        
        to_remove: list[tuple[str, str]] = []
        for u, v, data in list(h.edges(data=True)):
            bonds = data.get("bonds", [])
            kept = []
            for link in bonds:
                pair = frozenset((link.a1_symbol, link.a2_symbol))
                if pair in allow_pairs:
                    kept.append(link)

            if not kept:
                to_remove.append((u, v))
            else:
                data["bonds"] = kept
                data["n_bonds"] = len(kept)

        h.remove_edges_from(to_remove)

        if drop_isolated:
            iso = [n for n in h.nodes if h.degree(n) == 0 and n != self.unassigned]
            h.remove_nodes_from(iso)
        
        return AssemblyGraph(g=h, unassigned=self.unassigned, validate_upon_initialization=True)
    
    def connected_components(self, keep_unassigned: bool = False) -> list["AssemblyGraph"]:
        """
        Get the connected components of the assembly graph.

        :param keep_unassigned: whether to keep the unassigned node in components (default: False)
        :return: list of AssemblyGraph instances representing connected components
        """
        h = self.g.copy()

        if not keep_unassigned and h.has_node(self.unassigned):
            h.remove_node(self.unassigned)

        comps: list["AssemblyGraph"] = []
        for nodes in nx.connected_components(h):
            sub = h.subgraph(nodes).copy()
            comps.append(AssemblyGraph(g=sub, unassigned=self.unassigned, validate_upon_initialization=True))
        
        return comps
    
    def longest_path(self, keep_unassigned: bool = False, max_starts: int = 25) -> list[MolNode]:
        """
        Find the longest path of monomer nodes in assembly graph.

        :param keep_unassigned: whether to keep the unassigned node in the graph (default: False)
        :param max_starts: maximum number of starting nodes for greedy search (default: 25)
        :return: list of MolNode instances representing the longest path
        """
        g = self.g
        h = g.copy()

        if not keep_unassigned and hasattr(self, "unassigned") and h.has_node(self.unassigned):
            h.remove_node(self.unassigned)

        # Double check that every node has a MolNode attached
        for n in h.nodes:
            if "molnode" not in h.nodes[n] or h.nodes[n]["molnode"] is None:
                raise ValueError(f"AssemblyGraph node {n!r} has no valid 'molnode' attached")

        # Empty graph case
        if h.number_of_nodes() == 0:
            return []

        def node_to_molnode(node_id: Any) -> MolNode:
            """
            Convert a graph node ID to its corresponding MolNode.
            
            :param node_id: node ID in the graph
            :return: corresponding MolNode
            """
            mn = h.nodes[node_id].get("molnode", None)
            if mn is None:
                raise ValueError(f"Node {node_id!r} has no 'molnode' attached")

            return mn

        def to_molnodes(path_nodes: list[Any]) -> list[MolNode]:
            """
            Convert a list of graph node IDs to their corresponding MolNodes.
            
            :param path_nodes: list of node IDs in the graph
            :return: list of corresponding MolNodes
            """
            return [node_to_molnode(n) for n in path_nodes]


        for comp_nodes in nx.connected_components(h):
            # Work per connected component; pick the longest result
            best_path_nodes: list[Any] = []
            
            c = h.subgraph(comp_nodes).copy()
            if c.number_of_nodes() == 0:
                continue

            # Monomer graph is a tree: use diameter via two BFS
            is_tree = nx.is_connected(c) and (c.number_of_edges() == c.number_of_nodes() - 1)
            if is_tree:
                # Diameter of a tree via two BFS
                start = next(iter(c.nodes))
                dist1 = nx.single_source_shortest_path_length(c, start)
                u = max(dist1, key=dist1.get)

                dist2 = nx.single_source_shortest_path_length(c, u)
                v = max(dist2, key=dist2.get)

                path_uv = nx.shortest_path(c, u, v)
                if len(path_uv) > len(best_path_nodes):
                    best_path_nodes = path_uv
                continue

            # General/cyclic case: multi-start + 1-step lookahead
            def greedy_walk(start: Any) -> list[Any]:
                """
                Perform a greedy walk starting from the given node.
                
                :param start: starting node ID
                :return: list of node IDs in the greedy path
                """
                used = {start}
                path = [start]
                cur = start

                while True:
                    options = [nb for nb in c.neighbors(cur) if nb not in used]
                    if not options:
                        break

                    best_nb = None
                    best_score = None

                    # Precompute allowed nodes list once per step
                    # (unvisited nodes plus candidate neighbor)
                    all_nodes = list(c.nodes)

                    for nb in options:
                        allowed = [x for x in all_nodes if x not in used] + [nb]
                        sub = c.subgraph(allowed)

                        # How many nodes remain reachable if we go to nb?
                        reachable = nx.single_source_shortest_path_length(sub, nb)
                        score = (len(reachable), -c.degree(nb))  # prefer more reachable, tie-break lower degree

                        if best_score is None or score > best_score:
                            best_score = score
                            best_nb = nb

                    if best_nb is None:
                        break

                    used.add(best_nb)
                    path.append(best_nb)
                    cur = best_nb

                return path

            # Try low-degree starts first (often good for long simple paths)
            nodes_sorted = sorted(c.nodes, key=lambda n: c.degree(n))
            n_starts = min(max_starts, len(nodes_sorted))
            starts = nodes_sorted[:n_starts]

            best_comp_path: list[Any] = []
            for s in starts:
                cand = greedy_walk(s)
                if len(cand) > len(best_comp_path):
                    best_comp_path = cand

            # Try starting from the found endpoints (sometimes extends)
            if best_comp_path:
                for endpoint in (best_comp_path[0], best_comp_path[-1]):
                    cand = greedy_walk(endpoint)
                    if len(cand) > len(best_comp_path):
                        best_comp_path = cand

            if len(best_comp_path) > len(best_path_nodes):
                best_path_nodes = best_comp_path

        return to_molnodes(best_path_nodes)
    

    def iter_kmers(
        self,
        k: int,
        include_unassigned: bool = False,
        identified_only: bool = False
    ) -> Generator[tuple[MolNode, ...], None, None]:
        """
        Iterate over all k-mers (paths of length k) in the assembly graph.

        :param k: length of the k-mers (must be at least 1)
        :param include_unassigned: whether to include the unassigned node in paths (default: False)
        :param identified_only: whether to yield only k-mers with all identified monomers (default: False)
        :yield: tuples of MolNode instances representing k-mers
        """
        if k < 1:
            raise ValueError("k must be at least 1")
        
        g = self.g

        def usable_node_ids() -> list[str]:
            ids = []
            for n in g.nodes:
                if (not include_unassigned) and (n == self.unassigned):
                    continue
                ids.append(n)
            return ids
        
        def node_to_molnode(node_id: str) -> MolNode:
            mn = g.nodes[node_id].get("molnode", None)
            if mn is None:
                raise ValueError(f"AssemblyGraph node {node_id!r} has no 'molnode' attached")
            return mn
        
        node_ids = usable_node_ids()

        # k == 1: one k-mer per node
        if k == 1:
            for nid in node_ids:
                mn = node_to_molnode(nid)
                if identified_only and not mn.is_identified:
                    continue
                yield (mn,)
            return
        
        # Stack items are (current_node_id, path_node_ids)
        stack: list[tuple[str, list[str]]] = [(start, [start]) for start in node_ids]

        while stack:
            cur, path = stack.pop()

            if len(path) == k:
                km = tuple(node_to_molnode(pid) for pid in path)
                if identified_only and any(not n.is_identified for n in km):
                    continue
                yield km
                continue

            for nbr in g.neighbors(cur):
                if (not include_unassigned) and (nbr == self.unassigned):
                    continue
                stack.append((nbr, path + [nbr]))

    
    def validate(self) -> None:
        """
        Validate the assembly graph structure.
        """
        for n, data in self.g.nodes(data=True):
            for k in ("molnode", "tags", "identity"):
                if k not in data:
                    raise ValueError(f"AssemblyGraph node {n!r} missing required attribute {k!r}")
                
            if not isinstance(data["tags"], set):
                raise ValueError(f"AssemblyGraph node {n!r} tags must be set[int]")
            
            if data["molnode"] is None and n != self.unassigned:
                raise ValueError(f"AssemblyGraph node {n!r} has None molnode but is not unassigned node")
            
        for u, v, data in self.g.edges(data=True):
            for k in ("bonds", "n_bonds"):
                if k not in data:
                    raise ValueError(f"AssemblyGraph edge {u!r}-{v!r} missing required attribute {k!r}")
                
            if not isinstance(data["bonds"], list):
                raise ValueError(f"AssemblyGraph edge {u!r}-{v!r} bonds must be list[RootBondLink]")
            
            if not isinstance(data["n_bonds"], int):
                raise ValueError(f"AssemblyGraph edge {u!r}-{v!r} n_bonds must be int")
            
    def to_dict(self) -> dict[str, Any]:
        """
        Convert the AssemblyGraph to a dictionary.

        :return: dictionary representation of the AssemblyGraph
        """
        nodes_out: list[dict[str, Any]] = []
        for node_id, data in self.g.nodes(data=True):
            tags = data.get("tags", set())
            tags_json = sorted(tags)  # stable + JSON-friendly

            mn = data.get("molnode", None)
            mn_json = None if mn is None else mn.to_dict()

            ident = data.get("identity", None)
            ident_json = None if ident is None else ident.to_dict()

            nodes_out.append(
                {
                    "id": node_id,
                    "tags": tags_json,
                    "identity": ident_json,
                    "molnode": mn_json,
                }
            )

        edges_out: list[dict[str, Any]] = []
        for u, v, data in self.g.edges(data=True):
            bonds = data.get("bonds", [])
            edges_out.append(
                {
                    "u": u,
                    "v": v,
                    "bonds": [b.to_dict() for b in bonds],
                    "n_bonds": int(data.get("n_bonds", len(bonds))),
                }
            )

        return {
            "unassigned": self.unassigned,
            "nodes": nodes_out,
            "edges": edges_out,
        }
    
    @classmethod
    def from_dict(cls, data: dict[str, Any], validate: bool = True) -> "AssemblyGraph":
        """
        Create an AssemblyGraph from a dictionary.

        :param data: dictionary representation of the AssemblyGraph
        :param validate: whether to validate the graph after creation (default: True)
        :return: AssemblyGraph instance
        """
        unassigned = data.get("unassigned", "unassigned")
        g = nx.Graph()

        # Nodes
        for nd in data.get("nodes", []):
            node_id = nd["id"]
            tags = set(nd.get("tags", []))

            mn_payload = nd.get("molnode", None)
            molnode = None if mn_payload is None else MolNode.from_dict(mn_payload)

            ident_payload = nd.get("identity", None)
            identity = None if ident_payload is None else MolIdentity.from_dict(ident_payload)

            g.add_node(node_id, molnode=molnode, tags=tags, identity=identity)

        # Edges
        for ed in data.get("edges", []):
            u = ed["u"]
            v = ed["v"]

            bonds_raw = ed.get("bonds", [])
            bonds = [RootBondLink.from_dict(b) for b in bonds_raw]

            n_bonds = int(ed.get("n_bonds", len(bonds)))
            g.add_edge(u, v, bonds=bonds, n_bonds=n_bonds)

        ag = cls(g=g, unassigned=unassigned, validate_upon_initialization=False)

        if validate:
            ag.validate()
            
        return ag

    @classmethod
    def build(
        cls,
        root_mol: Mol,
        monomers: Iterable[MolNode],
        allow_overlaps: bool = False,
        include_unassigned: bool = False,
        unassigned: str = "unassigned",
        validate: bool = True,
    ) -> "AssemblyGraph":
        """
        Build an AssemblyGraph from the given root molecule and monomers.

        :param root_mol: RDKit Mol representing the root molecule
        :param monomers: iterable of MolNode representing the monomers
        :param allow_overlaps: whether to allow overlapping monomers (default: False)
        :param include_unassigned: whether to include unassigned regions as a node (default: False)
        :param unassigned: name of the unassigned node (default: "unassigned")
        :param validate: whether to validate the graph after building (default: True)
        :return: AssemblyGraph instance
        """
        g = build_assembly_graph(
            root_mol=root_mol,
            monomers=monomers,
            allow_overlaps=allow_overlaps,
            include_unassigned=include_unassigned,
        )
        ag = cls(g=g, unassigned=unassigned)

        if validate:
            ag.validate()

        return ag
    
    def draw(
        self,
        with_labels: bool = True,
        show_unassigned: bool = False,
        node_size: int = 1600,
        font_size: int = 9,
        edge_with_scale: float = 1.0,
        savepath: str | None = None,
    ) -> None:
        """
        Visualize the assembly graph using Matplotlib.

        :param with_labels: whether to show node labels (default: True)
        :param show_unassinged: whether to show the unassigned node (default: False)
        :param node_size: size of the nodes (default: 1600)
        :param font_size: font size for labels (default: 9)
        :param edge_with_scale: scale factor for edge widths (default: 1.0)
        :param savepath: optional path to save the figure (default: None)
        """
        # Hide unassigned if requested
        if not show_unassigned:
            g = self.drop_unassigned().g
        else:
            g = self.g

        if g.number_of_nodes() == 0:
            raise ValueError("AssemblyGraph has no nodes to show")
        
        # Layout
        pos = nx.spring_layout(g, seed=42)

        # Node colors
        node_colors = []
        labels = {}

        for n, data in g.nodes(data=True):
            if n == self.unassigned:
                node_colors.append("lightgray")
                labels[n] = "unassigned"
                continue

            ident = data.get("identity")
            if ident is None:
                node_colors.append("lightblue")
                labels[n] = n[:8]
            else:
                node_colors.append("lightgreen")
                labels[n] = getattr(ident, "name", n[:8])

        # Edge widths from number of root bonds
        widths = [edge_with_scale * max(1, data.get("n_bonds", 1)) for _, _, data in g.edges(data=True)]

        plt.figure(figsize=(8, 8))
        nx.draw_networkx_nodes(
            g,
            pos,
            node_color=node_colors,
            node_size=node_size,
            edgecolors="black",
        )
        nx.draw_networkx_edges(
            g,
            pos,
            width=widths,
            alpha=0.8,
        )
        
        if with_labels:
            nx.draw_networkx_labels(
                g,
                pos,
                labels=labels,
                font_size=font_size,
            )

        plt.axis("off")
        plt.tight_layout()

        if savepath is not None:
            plt.savefig(savepath, dpi=300)
        else:
            plt.show()
