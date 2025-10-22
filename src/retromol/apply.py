"""This module contains functions for applying custom rules to molecules."""

import itertools
import logging
from copy import deepcopy
from dataclasses import dataclass
from typing import Any

from networkx import Graph

from retromol import matching, rules
from retromol.chem import (
    ChemicalReaction,
    Mol,
    encode_mol,
    get_tags_mol,
    mol_to_smiles,
    neutralize_mol,
    smiles_to_mol,
)
from retromol.config import LOGGER_NAME
from retromol.errors import MotifGraphNodeWithoutAttributesError
from retromol.graph import merge_nodes, mol_to_graph
from retromol.io import Input as RetroMolInput
from retromol.rules import DummyReactionRule, ReactionRule


def _reactive_template_atoms(rxn: ChemicalReaction) -> list[set[int]]:
    """
    For each reactant-template in rxn, return the set of template-atom-indices
    that actually change (i.e. have a broken/formed bond or disappear/appear).
    We return a list: one set per reactant-template in the order they appear.

    :param rxn: RDKit ChemicalReaction object
    :return: List of sets, each set contains indices of reactive atoms in the corresponding reactant template
    """
    # First, build a map from map‐no -> (reactant_template_idx, reactant_atom_idx)
    reactant_maps: dict[
        int, tuple[int, int]
    ] = {}  # map_no -> (which reactant‐template, which atom‐idx in that template)
    for ri in range(rxn.GetNumReactantTemplates()):
        templ = rxn.GetReactantTemplate(ri)
        for atom in templ.GetAtoms():
            mnum = atom.GetAtomMapNum()
            if mnum:
                reactant_maps[mnum] = (ri, atom.GetIdx())

    # Next, build a map from map‐no -> (which product_template_idx, product_atom_idx)
    product_maps: dict[int, tuple[int, int]] = {}
    for pi in range(rxn.GetNumProductTemplates()):
        templ_p = rxn.GetProductTemplate(pi)
        for atom in templ_p.GetAtoms():
            mnum = atom.GetAtomMapNum()
            if mnum:
                product_maps[mnum] = (pi, atom.GetIdx())

    # Now we scan each reactant‐template atom and see if it "persists" into product with the same adjacency,
    # or if its bonding pattern changes, or if it disappears entirely. If any of those are true -> it's reactive.
    reactive_sets: list[set[int]] = [set() for _ in range(rxn.GetNumReactantTemplates())]

    # Pre‐compute adjacency‐lists (by map‐number) for reactant vs. product
    #  – build map_no -> set(of neighbor‐map_numbers) in reactant and product
    react_adj: dict[int, set[int]] = {}
    prod_adj: dict[int, set[int]] = {}

    # Build reactant adjacency by map‐num
    for ri in range(rxn.GetNumReactantTemplates()):
        templ = rxn.GetReactantTemplate(ri)
        for bond in templ.GetBonds():
            a1, a2 = bond.GetBeginAtom(), bond.GetEndAtom()
            m1, m2 = a1.GetAtomMapNum(), a2.GetAtomMapNum()
            if m1 and m2:
                react_adj.setdefault(m1, set()).add(m2)
                react_adj.setdefault(m2, set()).add(m1)

    # Build product adjacency by map‐num
    for pi in range(rxn.GetNumProductTemplates()):
        templ_p = rxn.GetProductTemplate(pi)
        for bond in templ_p.GetBonds():
            a1_p, a2_p = bond.GetBeginAtom(), bond.GetEndAtom()
            m1, m2 = a1_p.GetAtomMapNum(), a2_p.GetAtomMapNum()
            if m1 and m2:
                prod_adj.setdefault(m1, set()).add(m2)
                prod_adj.setdefault(m2, set()).add(m1)

    # Now: for each map_no in the reactant_templates, check:
    #  (a) if that map_no does NOT appear in product_maps at all -> the atom was deleted (= reactive)
    #  (b) if it DOES appear, compare react_adj[map_no] vs. prod_adj[map_no].
    #      If they differ -> bond‐pattern changed -> reactive
    #  (c) also check if atomic number or formal charge changed (rare in a template, but could).
    #      We compare the two atoms directly. We need to find the reactant‐template Atom and product‐template
    #      Atom to compare.
    for mnum, (rtempl_idx, ratom_idx) in reactant_maps.items():
        if mnum not in product_maps:
            # Disappeared in the product – this atom is definitely reactive
            reactive_sets[rtempl_idx].add(ratom_idx)
        else:
            # Compare adjacency
            react_neighbors = react_adj.get(mnum, set())
            prod_neighbors = prod_adj.get(mnum, set())
            if react_neighbors != prod_neighbors:
                reactive_sets[rtempl_idx].add(ratom_idx)
            else:
                # Check if element or charge changed
                (pi, patom_idx) = product_maps[mnum]
                react_atom = rxn.GetReactantTemplate(rtempl_idx).GetAtomWithIdx(ratom_idx)
                prod_atom = rxn.GetProductTemplate(pi).GetAtomWithIdx(patom_idx)
                if (
                    react_atom.GetAtomicNum() != prod_atom.GetAtomicNum()
                    or react_atom.GetFormalCharge() != prod_atom.GetFormalCharge()
                ):
                    # If neither bonding‐pattern nor element‐/charge changed, it is NOT reactive
                    reactive_sets[rtempl_idx].add(ratom_idx)

    return reactive_sets


def _index_uncontested(
    mol: Mol,
    rls: list[ReactionRule],
    failed_combos: set[tuple[int, frozenset[int]]],
) -> list[tuple[ReactionRule, set[int]]]:
    """
    Index uncontested reactions for applying preprocessing rules in bulk.

    :param mol: RDKit molecule
    :param rls: List of preprocessing rules
    :param failed_combos: Set of failed combinations to avoid infinite loops
    :return: Uncontested reactions
    """
    up_for_election: list[tuple[ReactionRule, set[int], set[int]]] = []
    for rl in rls:
        if not rl.rxn:
            continue  # skip rules without a reaction template

        reactive_inds = _reactive_template_atoms(rl.rxn)[0]
        all_reactant_matches: list[tuple[tuple[int, ...], ...]] = []
        all_reactant_matches_reactive_items: list[list[list[int]]] = []
        for template_idx in range(rl.rxn.GetNumReactantTemplates()):
            reactant_template = rl.rxn.GetReactantTemplate(template_idx)
            reactant_matches: tuple[tuple[int, ...], ...] = mol.GetSubstructMatches(reactant_template)
            all_reactant_matches.append(reactant_matches)
            new_reactant_matches: list[list[int]] = []
            for reactant_match in reactant_matches:
                new_reactant_matches.append([reactant_match[idx] for idx in reactive_inds])
            all_reactant_matches_reactive_items.append(new_reactant_matches)

        # Generate all possible match sets, for when reaction template matches multiple sites
        match_sets = list(itertools.product(*all_reactant_matches))
        match_sets_reactive_items = list(itertools.product(*all_reactant_matches_reactive_items))
        match_sets = [set(itertools.chain(*match_set)) for match_set in match_sets]
        match_sets_reactive_items = [set(itertools.chain(*match_set)) for match_set in match_sets_reactive_items]
        for match_set, match_set_reactive_items in zip(match_sets, match_sets_reactive_items, strict=True):
            up_for_election.append((rl, match_set, match_set_reactive_items))

    # Check which reactions with matched templates are uncontested and which are contested
    uncontested: list[tuple[ReactionRule, set[int]]] = []
    for i, (rl, match_set, match_set_reactive_items) in enumerate(up_for_election):
        # Rules with ring matching conditions are always contested
        if rl.has_ring_matching_condition():
            continue

        # Check if match set has overlap with any other match set
        # has_overlap = any(match_set.intersection(o) for j, (_, o, o_r) in enumerate(up_for_election) if i != j)
        has_overlap = any(
            match_set_reactive_items.intersection(o_r) for j, (_, _, o_r) in enumerate(up_for_election) if i != j
        )
        if not has_overlap:
            uncontested.append((rl, match_set))

    # Filter out failed combinations to avoid infinite loops
    uncontested = [
        (rl, match_set) for rl, match_set in uncontested if (rl.id, frozenset(match_set)) not in failed_combos
    ]

    return uncontested


def _apply_uncontested(
    parent: Mol,
    uncontested: list[tuple[rules.ReactionRule, set[int]]],
    original_taken_tags: list[int],
) -> tuple[list[Mol], set[str], set[tuple[int, frozenset[int]]]]:
    """
    Apply uncontested reactions in bulk.

    :param parent: RDKit molecule
    :param uncontested: List of uncontested reactions
    :param original_taken_tags: List of atom tags from original reactant
    :return: List of trtue products, a set of applied reaction ids,  and a set of failed combinations
    """
    applied_reactions: set[str] = set()

    tags_in_parent: set[int] = set(get_tags_mol(parent))

    # We make sure all atoms, even the ones not from original reactant, have a
    # unique isotope number, so we can track them through consecutive reactions
    temp_taken_tags = get_tags_mol(parent)
    for atom in parent.GetAtoms():
        if atom.GetIsotope() == 0:
            tag = 1
            while tag in original_taken_tags or tag in temp_taken_tags:
                tag += 1
            atom.SetIsotope(tag)
            temp_taken_tags.append(tag)

    # Validate that all atoms have a unique tag
    num_tagged_atoms = len(set(get_tags_mol(parent)))
    if num_tagged_atoms != len(parent.GetAtoms()):
        raise ValueError("Not all atoms have a unique tag before applying uncontested reactions")

    # Map tags to atomic nums so we can create masks and reassign atomic nums later on
    idx_to_tag = {a.GetIdx(): a.GetIsotope() for a in parent.GetAtoms()}

    # All uncontested reactions become a single node in the reaction_graph
    products: list[Mol] = []
    failed_combos: set[tuple[int, frozenset[int]]] = set()  # keep track of failed combinations to avoid infinite loops

    for rl, match_set in uncontested:
        msk = set([idx_to_tag[idx] for idx in match_set])  # create mask for reaction

        # We use the input parent if there are no products, otherwise we have to find out
        # which product now contains the mask (i.e., the reaction template for this reaction)
        if len(products) != 0:
            new_parent: Mol | None = None
            for product in products:
                product_tags = set(get_tags_mol(product))
                if msk.issubset(product_tags):
                    new_parent = product
                    products = [p for p in products if p != product]
                    break

            if new_parent is None:
                # raise ValueError("no product found that contains the mask")
                # If no product is found, we continue with the next uncontested reaction
                continue

            parent = new_parent

        # Register all tags currently taken by atoms in parent
        temp_taken_tags_uncontested = get_tags_mol(parent)

        # Newly introduced atoms by one of the uncontested reactions need a unique tag
        for atom in parent.GetAtoms():
            if atom.GetIsotope() == 0:  # newly introduced atom has tag 0
                # Loop until we find a tag that is not already taken
                tag = 1
                while tag in (temp_taken_tags_uncontested + original_taken_tags + temp_taken_tags):
                    tag += 1
                atom.SetIsotope(tag)
                temp_taken_tags_uncontested.append(tag)

        unmasked_parent = deepcopy(parent)  # keep original parent for later
        results = rl(parent, msk)  # apply reaction rule

        try:
            if len(results) == 0:
                raise ValueError(f"No products from uncontested reaction {rl.rid}")

            if len(results) > 1:
                raise ValueError(f"More than one product from uncontested reaction {rl.rid}")

            result = results[0]
            applied_reactions.add(rl.rid)  # keep track of successfully applied reactions

            # Reset atom tags in products for atoms not in original reactant
            for product in result:
                for atom in product.GetAtoms():
                    if atom.GetIsotope() not in original_taken_tags and atom.GetIsotope() != 0:
                        atom.SetIsotope(0)
                products.append(product)

        except Exception:
            # Start function again with the next uncontested reaction
            for atom in parent.GetAtoms():
                if atom.GetIsotope() not in original_taken_tags and atom.GetIsotope() != 0:
                    atom.SetIsotope(0)
            products.append(unmasked_parent)
            failed_combos.add(
                (
                    rl.id,
                    frozenset(match_set),
                )
            )

    for product in products:
        # Any tag in product that is not in parent should be 0; otherwise we run into issues with
        # the set cover algorithm
        for atom in product.GetAtoms():
            if atom.GetIsotope() not in tags_in_parent and atom.GetIsotope() != 0:
                atom.SetIsotope(0)

    return products, applied_reactions, failed_combos


@dataclass
class ProcessingResult:
    """
    Data structure for processed data.

    :param enc_to_mol: Maps mol hash to Chem.Mol
    :param enc_to_rxn: Maps rxn index to Reaction
    :param rxn_graph: As {reactant_mol_encoding: {reaction_encoding: [child_mol_encodings], ...}, ...}
    :param applied_rxns: Set of successfully applied reaction ids
    """

    enc_to_mol: dict[str, Mol]  # encoding is a canonical SMILES with isotopic tags
    enc_to_rxn: dict[int, rules.ReactionRule]
    rxn_graph: dict[str, dict[int, list[str]]]
    applied_rxns: set[str]


class Input:
    """
    This class describes the input for a RetroMol run.
    """

    def __init__(self, cid: str, smi: str) -> None:
        """
        Initialize the input compound.

        :param cid: Compound identifier
        :param smi: SMILES representation
        """
        self.cid = cid

        # convert SMILES into RDKit molecule
        self.mol = smiles_to_mol(smi)
        neutralize_mol(self.mol)

        # store original atom indices as isotope number
        for atom in self.mol.GetAtoms():
            idx = atom.GetIdx()
            tag = idx + 1
            atom.SetIsotope(tag)

        # store SMILES representation with tags
        self._smi = mol_to_smiles(self.mol)

    def get_tags(self) -> list[int]:
        """
        Get the atom tags.

        :return: atom tags
        """
        return get_tags_mol(self.mol)


def process_mol(inp: RetroMolInput, reaction_rules: list[rules.ReactionRule]) -> ProcessingResult:
    """
    Apply custom rules to linearize a SMILES string.

    :param inp: Input object
    :param reaction_rules: List of processing rules
    :param reserved_nodes: Set of atom tags that should not be used for remapping
    :return: Preprocessed data structures
    """
    logger = logging.getLogger(LOGGER_NAME)

    # Setup processing
    original_taken_tags = inp.get_tags()
    rst_pre = ProcessingResult({}, {}, {}, set())
    num_rxn_nodes = 0  # keeps track of number of reaction nodes in reaction graph
    mols = [deepcopy(inp.mol)]  # queue of molecules to process

    # Set of (rule id, frozenset of match_set) tuples to track failed combinations
    failed_combos: set[tuple[int, frozenset[int]]] = set()  # keep track of failed combinations to avoid infinite loops

    # Process queue
    while mols:
        parent = mols.pop(0)

        # Encode parent molecule
        parent_encoding = encode_mol(parent)
        if parent_encoding not in rst_pre.enc_to_mol:
            rst_pre.enc_to_mol[parent_encoding] = deepcopy(parent)
            rst_pre.rxn_graph[parent_encoding] = {}

        # Index uncontested reactions
        uncontested = _index_uncontested(parent, reaction_rules, failed_combos)

        # Apply uncontested reactions in bulk
        if uncontested:
            products, applied_in_bulk, new_failed_combos = _apply_uncontested(parent, uncontested, original_taken_tags)
            logger.debug(
                f"Uncontested reactions applied ({len(applied_in_bulk)}): {', '.join([rl.rid for rl, _ in uncontested])}"  # noqa E501
            )
            logger.debug(
                f"Uncontested reactions failed ({len(new_failed_combos)}): {', '.join([str(id) for id, _ in new_failed_combos])}"  # noqa E501
            )
            failed_combos.update(new_failed_combos)

            # If all uncontested reactions failed, we continue with the next parent
            logger.debug(f"Found {len(products)} product(s) from uncontested reactions")

            # If no reaction was successful, we continue with the next parent
            if len(applied_in_bulk) == 0:
                logger.debug("No uncontested reactions applied, continuing with next parent")
                # Apparently everything failed, so we need to reparse the parent and try again with contested reactions
                mols.append(parent)
                continue

            if products:
                rst_pre.applied_rxns.update(applied_in_bulk)

                # All products are now products of our combined reaction node that contains all uncontested reactions
                num_rxn_nodes += 1
                rst_pre.enc_to_rxn[num_rxn_nodes] = DummyReactionRule("uncontested")
                rst_pre.rxn_graph[parent_encoding][num_rxn_nodes] = list()
                for product in products:
                    product_encoding = encode_mol(product)
                    if product_encoding not in rst_pre.rxn_graph[parent_encoding][num_rxn_nodes]:
                        rst_pre.rxn_graph[parent_encoding][num_rxn_nodes].append(product_encoding)
                    if product_encoding not in rst_pre.enc_to_mol:
                        rst_pre.enc_to_mol[product_encoding] = deepcopy(product)
                        rst_pre.rxn_graph[product_encoding] = {}
                        mols.append(product)

                # Restart loop with new parent
                continue

        # Exhaustive reaction_rule application for all contested reactions
        for rl in reaction_rules:
            results = rl(parent, None)  # apply reaction rule

            if results:
                logger.debug(f"Contested reaction {rl.rid} applied")

            for result in results:
                rst_pre.applied_rxns.add(rl.rid)  # keep track of successfully applied reactions

                # Encode reaction node
                num_rxn_nodes += 1
                if num_rxn_nodes not in rst_pre.enc_to_rxn:
                    rst_pre.enc_to_rxn[num_rxn_nodes] = rl
                else:
                    raise ValueError(f"reaction node {num_rxn_nodes} already exists for reaction {rl.rid}")

                # Encode product molecules
                rst_pre.rxn_graph[parent_encoding][num_rxn_nodes] = list()
                for child in result:
                    child_encoding = encode_mol(child)
                    if child_encoding not in rst_pre.rxn_graph[parent_encoding][num_rxn_nodes]:
                        rst_pre.rxn_graph[parent_encoding][num_rxn_nodes].append(child_encoding)

                    # Add child to queue if not already present in encoding_to_mol (i.e., previously processed)
                    if child_encoding not in rst_pre.enc_to_mol:
                        rst_pre.enc_to_mol[child_encoding] = deepcopy(child)
                        rst_pre.rxn_graph[child_encoding] = {}
                        mols.append(child)

    return rst_pre


def resolve_mol(
    mol: Mol,
    reserved_tags: set[int],
    reaction_rules: list[rules.ReactionRule],
    matching_rules: list[rules.MatchingRule],
    match_stereochemistry: bool,
    wave_config: dict[str, Any],
) -> "Graph[int | str]":
    """
    Apply custom rules to sequence a molecule into motif codes.

    :param mol: RDKit molecule
    :param reserved_tags: set of atom tags that should not be used for remapping
    :param reaction_rules: list of sequencing rules
    :param matching_rules: matching rules for identifying nodes
    :param match_stereochemistry: whether to match stereochemistry
    :param wave_config: configuration for the current wave
    :return: motif graph with merged and possibly identified nodes
    """
    logger = logging.getLogger(LOGGER_NAME)

    # Retrieve wave configuration
    use_leafs_only = wave_config.get("only_leaf_nodes", True)
    matching_groups = wave_config.get("matching_groups", None)

    # Filter matching rules based on priorities if provided
    if matching_groups is not None:
        matching_rules = [mr for mr in matching_rules if any([True for gr in mr.groups if gr in matching_groups])]
        logger.debug(f"Filtered matching rules based on priorities to {len(matching_rules)} rules")

    # Tag atoms without a tag yet with unique tags
    mol_to_process = RetroMolInput(
        "mol", mol, tag_compound=True, reserved_tags=deepcopy(reserved_tags)
    )  # copy because we modify it in place
    all_tags = mol_to_process.get_tags()
    all_tags.sort()

    logger.debug(f"Processing molecule {mol_to_process.cid} with SMILES {mol_to_process.smi}")
    logger.debug(f"Tags in molecule: {all_tags}")

    processing_result = process_mol(mol_to_process, reaction_rules)

    # Decide what nodes to use from processing result
    if use_leafs_only:
        # Only use leaf nodes from the reaction graph
        encoding_to_mol = {
            enc: mol
            for enc, mol in processing_result.enc_to_mol.items()
            if not processing_result.rxn_graph.get(enc)  # check if node has children
        }

        # If there are no leaf nodes, resort to all nodes
        # This might happen when an uncontested rule fails and returns itself as product
        if not encoding_to_mol:
            encoding_to_mol = processing_result.enc_to_mol

    else:
        # Use all nodes from the reaction graph
        encoding_to_mol = processing_result.enc_to_mol

    # Identify nodes and pick best set of identified nodes
    encoding_to_mol_identified = matching.identify_nodes(encoding_to_mol, matching_rules, match_stereochemistry)
    identified_nodes = list(encoding_to_mol_identified.keys())  # keys are the identified nodes

    unidentified_nodes = [node for node in encoding_to_mol.keys() if node not in identified_nodes]
    best_fit_identified, best_fit_unidentified = matching.solve_exact_cover_with_priority(
        encoding_to_mol, identified_nodes, unidentified_nodes, all_tags
    )

    # Create monomer graph for best fit for identified and unidentified nodes
    motif_graph: Graph[int | str] = mol_to_graph(mol_to_process.mol, use_tags=True)

    for node in best_fit_identified + best_fit_unidentified:
        # Merged nodes have a str key, unmerged nodes have an int key
        if not isinstance(node, str):
            continue  # skip unmerged nodes
        node_mol = encoding_to_mol[node]
        node_smiles = mol_to_smiles(node_mol)
        tags = get_tags_mol(node_mol)
        merge_nodes(motif_graph, merged_node_id=node, nodes=tags)  # modifies node, merged node has str key

        # Get identity of node
        if identity_with_props := encoding_to_mol_identified.get(node, None):
            identity = identity_with_props["identity"]
            props = identity_with_props["props"]
        else:
            identity = None
            props = {}

        logger.debug(f"Node {node} ({node_smiles}) has identity {identity}")

        # Make sure node_smiles only has reserved tags, if any, for mapping in readout functions
        node_mol = smiles_to_mol(node_smiles)
        for atom in node_mol.GetAtoms():
            tag = atom.GetIsotope()
            if tag not in reserved_tags:
                atom.SetIsotope(0)
        node_smiles = mol_to_smiles(node_mol)
        tags = get_tags_mol(node_mol)

        # Give properties to monomer graph nodes
        motif_graph.nodes[node]["graph"] = None
        motif_graph.nodes[node]["identity"] = identity
        motif_graph.nodes[node]["props"] = props
        motif_graph.nodes[node]["tags"] = get_tags_mol(node_mol)
        motif_graph.nodes[node]["smiles"] = node_smiles
        motif_graph.nodes[node]["smiles_no_tags"] = mol_to_smiles(node_mol, remove_tags=True)
        motif_graph.nodes[node]["wave_name"] = wave_config.get("wave_name", None)

    # Check if there are ny nodes that have no attrs (something went wrong...)
    nodes_without_attrs = [n for n, d in motif_graph.nodes(data=True) if not d]
    if nodes_without_attrs:
        # This happened sometimes before we explicitly checked for overlapping atom
        # tags (see rules.OverlappingAtomTagsError). Might be caused by some other type
        # of behavior from RDKit, so we explicitly check for it now.
        raise MotifGraphNodeWithoutAttributesError(f"Nodes without attributes found: {nodes_without_attrs}")

    return motif_graph
