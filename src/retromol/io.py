"""This module describes the basic output class for RetroMol."""

import json
from collections import Counter, defaultdict
from dataclasses import dataclass
from typing import Any

from networkx import Graph, node_link_data, node_link_graph

from retromol.chem import (
    Mol,
    get_tags_mol,
    mol_to_inchikey,
    mol_to_smiles,
    neutralize_mol,
    smiles_to_mol,
)


class Input:
    """
    This class describes the input for a RetroMol run.
    """

    def __init__(
        self,
        cid: str,
        repr: Mol | str,
        props: dict[str, Any] | None = None,
        tag_compound: bool = True,
        reserved_tags: set[int] | None = None,
    ) -> None:
        """
        Initialize the input compound.

        :param cid: compound identifier
        :param mol: RDKit molecule or SMILES string
        :param props: additional properties
        :param tag_compound: whether to tag the compound's atoms, existing tags will be preserved
        """
        self.cid = cid
        self.props = props

        if isinstance(repr, Mol):
            self.mol = repr
            if tag_compound:
                reserved = reserved_tags or set()
                for atom in self.mol.GetAtoms():
                    tag = atom.GetIsotope()
                    if tag in reserved:
                        continue
                    if tag == 0:
                        idx = atom.GetIdx()
                        tag = idx + 1
                        while tag in reserved:
                            tag += 1
                        atom.SetIsotope(tag)
                        reserved.add(tag)

            # Store SMILES representation with tags
            self.smi = mol_to_smiles(self.mol)

        else:
            smi = repr

            # Sanitize SMILES
            smi = smi.replace("[N]", "N")  # avoid parsing issues with RDKit

            # Convert SMILES into RDKit molecule
            self.mol = smiles_to_mol(smi, retain_largest_fragment=True)
            neutralize_mol(self.mol)

            # Store original atom indices as isotope number
            if tag_compound:
                reserved = reserved_tags or set()
                for atom in self.mol.GetAtoms():
                    tag = atom.GetIsotope()
                    if tag in reserved:
                        continue
                    if tag == 0:
                        idx = atom.GetIdx()
                        tag = idx + 1
                        while tag in reserved:
                            tag += 1
                        atom.SetIsotope(tag)
                        reserved.add(tag)

            # Store SMILES representation with tags
            self.smi = mol_to_smiles(self.mol)

    def get_tags(self) -> list[int]:
        """
        Get the atom tags.

        :return: atom tags
        """
        return get_tags_mol(self.mol)


@dataclass
class Result:
    """
    A class representing the result of a RetroMol operation.

    :param input_id: a unique identifier for the input molecule
    :param graph: a networkx Graph representing the motif graph of the input molecule. Nodes of this graph may have
        a “graph” attribute which is itself another nx.Graph (or None), arbitrarily nested
    :param props: additional properties associated with the input molecule.
    :param sha256_reaction_rules: SHA256 hash of the reaction rules used (optional)
    :param sha256_matching_rules: SHA256 hash of the matching rules used (optional)
    """

    input_id: str
    graph: "Graph[int | str]"
    props: dict[str, Any] | None
    sha256_reaction_rules: str | None
    sha256_matching_rules: str | None

    def serialize(self) -> dict[str, Any]:
        """
        Serialize this Result to a JSON-friendly dict, including any nested graphs.

        :return: serialized representation of the Result
        """
        return {
            "input_id": self.input_id,
            "graph": self._serialize_graph(self.graph),
            "props": self.props,
            "sha256_reaction_rules": self.sha256_reaction_rules,
            "sha256_matching_rules": self.sha256_matching_rules,
        }

    def get_input_smiles(self, remove_tags: bool = False) -> str:
        """
        Get the SMILES representation of the input molecule.

        This is a convenience method to access the input's SMILES from the Result.

        :param remove_tags: if True, removes the atom tags from the SMILES.
        :return: SMILES string of the input molecule.
        """
        if remove_tags:
            smiles = self.graph.graph.get("smiles_no_tags", None)
        else:
            smiles = self.graph.graph.get("smiles", None)

        if smiles is None:
            raise ValueError("SMILES not found in graph attributes.")

        # Bit unsafe, but we trust our own data structure here
        return smiles

    def get_props(self) -> dict[str, Any]:
        """
        Get the additional properties associated with the input molecule.

        :return: dictionary of additional properties
        """
        return self.props if self.props is not None else {}

    @staticmethod
    def _serialize_graph(g: "Graph[int | str]") -> dict[str, Any]:
        # First use node_link_data to turn the graph structure into primitives + attrs
        data = node_link_data(g)

        # Now look for any node-attribute called "graph" that is itself an nx.Graph, and recurse
        for node in data["nodes"]:
            node["id"] = str(node["id"])  # Ensure node IDs are strings
            sub = node.get("graph")
            if isinstance(sub, Graph):
                node["graph"] = Result._serialize_graph(sub)
            # If it's None (or already a dict), leave it as is

        # Make sure to also stringify the link endpoints
        for link in data["links"]:
            link["source"] = str(link["source"])
            link["target"] = str(link["target"])

        return data

    @staticmethod
    def from_serialized(data: dict[str, Any]) -> "Result":
        """
        Reconstruct a Result from the dict form produced by serialize().

        :param data: serialized representation of the Result
        :return: reconstructed Result object
        """
        nested = data["graph"]
        g = Result._deserialize_graph(nested)
        props = data.get("props", None)
        sha256_reaction_rules = data.get("sha256_reaction_rules", None)
        sha256_matching_rules = data.get("sha256_matching_rules", None)
        return Result(
            input_id=data["input_id"],
            graph=g,
            props=props,
            sha256_reaction_rules=sha256_reaction_rules,
            sha256_matching_rules=sha256_matching_rules,
        )

    @staticmethod
    def _deserialize_graph(data: dict[str, Any]) -> "Graph[int | str]":
        # Convert string IDs back to integers
        for node in data["nodes"]:
            node["id"] = node["id"]

        # Convert string source/target IDs back to integers in links
        for link in data["links"]:
            link["source"] = link["source"]
            link["target"] = link["target"]

        # Build the Graph object from its node-link dict
        g = node_link_graph(data)

        # Walk each node, and if it has a "graph" attribute that is a dict, recurse
        for _, attrs in g.nodes(data=True):
            nested = attrs.get("graph")
            if isinstance(nested, dict):
                attrs["graph"] = Result._deserialize_graph(nested)

        return g

    def to_json(self) -> str:
        """
        Dump the fully serialized form to a JSON string.

        :return: JSON string representation of the Result
        """
        return json.dumps(self.serialize(), indent=2)

    @staticmethod
    def from_json(s: str) -> "Result":
        """
        Load a Result back from a JSON string.

        :param s: JSON string representation of the Result
        :return: reconstructed Result object
        """
        data = json.loads(s)
        return Result.from_serialized(data)

    def summarize_by_depth(self, propagate_through_identified: bool = True) -> dict[int, dict[str, Any]]:
        """
        For each depth d (root children are depth 1), accumulate UNIQUE input tags
        from all identified nodes seen at depths <= d.

        :param propagate_through_identified:
            if True, continue traversing into subgraphs of identified nodes
            if False, only traverse into unidentified nodes' subgraphs

        Returns { depth: {
            "covered_tag_count": int,      # cumulative unique tags
            "coverage": float,             # covered_tag_count / |input_tags|
            "n_nodes": int,                # per-depth
            "n_identified": int,           # per-depth
            "n_unidentified": int,         # per-depth
            "monomer_counts": Counter,     # per-depth (identity -> count)
            "wave_name": Optional[str],    # first seen at that depth
        }, ...}
        """
        # Input tag universe
        root_tags = self.graph.graph.get("tags")
        if isinstance(root_tags, (list, tuple, set)):
            T0 = {int(t) for t in root_tags}
        else:
            root_smi = self.graph.graph.get("smiles")
            if not root_smi:
                raise ValueError("Root graph missing 'tags' and 'smiles'.")
            T0 = set(get_tags_mol(smiles_to_mol(root_smi)))
        denom = len(T0) if T0 else 0

        # Per-depth accumulators
        tags_at_depth: dict[int, set[int]] = defaultdict(set)  # ONLY from identified nodes
        monomers_by_depth: dict[int, Counter[str]] = defaultdict(Counter)
        nodes_by_depth: dict[int, int] = defaultdict(int)
        identified_by_depth: dict[int, int] = defaultdict(int)
        unidentified_by_depth: dict[int, int] = defaultdict(int)
        wave_name_by_depth: dict[int, str] = {}

        def _walk(g, depth: int):
            nd = depth + 1
            for _, attrs in g.nodes(data=True):
                nodes_by_depth[nd] += 1

                identity = attrs.get("identity")
                node_tags_raw = attrs.get("tags")

                # Only identified nodes contribute tags
                if identity is not None:
                    identified_by_depth[nd] += 1
                    monomers_by_depth[nd][identity] += 1

                    if node_tags_raw is not None:
                        try:
                            node_tags = {int(t) for t in node_tags_raw}
                        except Exception:
                            node_tags = set()
                        # Only count tags that belong to the input
                        tags_at_depth[nd].update(node_tags & T0)
                else:
                    unidentified_by_depth[nd] += 1

                if nd not in wave_name_by_depth:
                    wn = attrs.get("wave_name")
                    if wn is not None:
                        wave_name_by_depth[nd] = wn

                sub = attrs.get("graph")
                if isinstance(sub, Graph):
                    if propagate_through_identified or identity is None:
                        _walk(sub, nd)

        _walk(self.graph, depth=0)

        # Cumulative union across depths
        depths: list[int] = sorted(set(nodes_by_depth) | set(tags_at_depth))
        cumulative: set[int] = set()
        summary: dict[int, dict[str, Any]] = {}
        for d in depths:
            cumulative |= tags_at_depth.get(d, set())
            covered = len(cumulative)
            cov = round((covered / denom), 4) if denom else 0.0
            summary[d] = {
                "wave_name": wave_name_by_depth.get(d, None),
                "monomer_counts": monomers_by_depth.get(d, Counter()),
                "covered_tag_count": covered,
                "coverage": cov,
                "n_nodes": nodes_by_depth.get(d, 0),
                "n_identified": identified_by_depth.get(d, 0),
                "n_unidentified": unidentified_by_depth.get(d, 0),
            }
        return summary

    def max_depth(self) -> int:
        """
        Return the deepest nesting level (wave index). Root is 0

        :return: Maximum depth of the nested motif graph
        """
        max_d = 0

        def _recurse(g: "Graph[int | str]", d: int):
            nonlocal max_d
            max_d = max(max_d, d)
            for _, attrs in g.nodes(data=True):
                sub = attrs.get("graph")
                if isinstance(sub, Graph):
                    _recurse(sub, d + 1)

        _recurse(self.graph, 0)
        return max_d

    def best_total_coverage(self, round_to: int = 2) -> float:
        """
        Best total coverage across all depths.

        :param round_to: number of decimal places to round the coverage to
        :return: best total coverage as a float rounded to 'round_to' decimal places
        """
        summary = self.summarize_by_depth(propagate_through_identified=False)
        total = max(info.get("coverage", 0.0) for info in summary.values())
        # clamp & round like calculate_coverage()
        total = max(0.0, min(1.0, total))
        return round(total, round_to)

    def get_unidentified_nodes(self) -> list[tuple[str, str]]:
        """
        Return all unidentified nodes as a list of (SMILES_without_tags, InChIKey).

        - Traverses the full nested motif graph
        - A node is considered unidentified if attrs.get("identity") is None
        - If a node already has "smiles_no_tags", uses it; otherwise removes
          isotopes/atom-maps from its "smiles" on the fly
        - Deduplicates and returns results sorted by SMILES, then InChIKey

        :return: list of tuples (SMILES_without_tags, InChIKey)
        """
        results: list[tuple[str, str]] = []
        seen: set[tuple[str, str]] = set()

        def _untag_smiles(tagged_smi: str) -> str:
            # Convert to RDKit mol, strip isotopes & atom map numbers, then re-SMILES
            mol = smiles_to_mol(tagged_smi)
            for atom in mol.GetAtoms():
                # Remove our tag carriers (we use isotopes for tags; atom-maps just in case)
                atom.SetIsotope(0)
                try:
                    atom.SetAtomMapNum(0)  # no-op if not set / older RDKit
                except Exception:
                    pass
            # Normalize charges
            neutralize_mol(mol)
            return mol_to_smiles(mol)

        def _inchikey_from_smiles(smi: str) -> str:
            """
            Get InChIKey from SMILES string.

            :param smi: SMILES string
            """
            mol = smiles_to_mol(smi)
            return mol_to_inchikey(mol)  # assumes this exists in retromol.chem

        def _walk(g: "Graph[int | str]") -> None:
            """
            Recursive walker to find unidentified nodes.

            :param g: current graph to walk
            """
            for _, attrs in g.nodes(data=True):
                # Collect if unidentified
                if attrs.get("identity") is None:
                    tagged = attrs.get("smiles")
                    if tagged is None:
                        continue
                    smi_no_tags = attrs.get("smiles_no_tags") or _untag_smiles(tagged)
                    try:
                        ik = _inchikey_from_smiles(smi_no_tags)
                    except Exception:
                        # Fall back to empty key if InChIKey generation fails
                        ik = ""
                    pair = (smi_no_tags, ik)
                    if pair not in seen:
                        seen.add(pair)
                        results.append(pair)

                # Recurse into subgraphs regardless of parent identification
                sub = attrs.get("graph")
                if isinstance(sub, Graph):
                    _walk(sub)

        _walk(self.graph)
        results.sort(key=lambda x: (x[0], x[1]))
        return results

    def get_identified_nodes(self) -> set[tuple[str, str, tuple[int, ...]]]:
        """
        Traverse the nested motif graph and return all identified nodes.

        :return: each entry is (identity, tags)
        """
        results: dict[tuple[str, tuple[int, ...]], str] = {}  # SMILES might vary for same identity+tags

        def _walk(g: "Graph[int | str]") -> None:
            for _, attrs in g.nodes(data=True):
                identity = attrs.get("identity")
                if identity is not None:
                    # Get untagged SMILES
                    smiles = attrs.get("smiles", "")
                    mol = smiles_to_mol(smiles)
                    smiles = mol_to_smiles(mol, remove_tags=True)

                    tags_list: list[int] = attrs.get("tags", [])
                    tags_list.sort()
                    tags: tuple[int, ...] = tuple(tags_list)
                    results[(identity, tags)] = smiles

                sub = attrs.get("graph")
                if isinstance(sub, Graph):
                    _walk(sub)

        _walk(self.graph)

        reformatted_results = {(key[0], smiles, key[1]) for (key, smiles) in results.items()}

        return reformatted_results
