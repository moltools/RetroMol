"""Module for RDKit reaction utilities."""

import itertools

from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdChemReactions import ChemicalReaction, ReactionFromSmarts

from retromol.model.rules import ReactionRule
from retromol.chem.mol import copy_mol
from retromol.chem.tagging import get_tags_mol


def smarts_to_reaction(smarts: str, use_smiles: bool = False) -> ChemicalReaction:
    """
    Converts a SMARTS string to an RDKit reaction.

    :param smarts: the SMARTS string to convert
    :param use_smiles: whether to interpret the SMARTS as SMILES
    :return: the RDKit reaction
    :raises ValueError: if the SMARTS pattern is invalid
    """
    rxn = ReactionFromSmarts(smarts, useSmiles=use_smiles)

    if rxn is None:
        raise ValueError(f"invalid reaction SMARTS: {smarts}")

    return rxn


def reactive_template_atoms(rxn: ChemicalReaction) -> list[set[int]]:
    """
    For each reactant-template in rxn, return the set of template-atom-indices
    that actually change (i.e. have a broken/formed bond or disappear/appear).
    We return a list: one set per reactant-template in the order they appear.

    :param rxn: RDKit ChemicalReaction object
    :return: List of sets, each set contains indices of reactive atoms in the corresponding reactant template
    """
    # First, build a map from map‐no -> (reactant_template_idx, reactant_atom_idx)
    reactant_maps: dict[int, tuple[int, int]] = {}  # map_no -> (which reactant‐template, which atom‐idx in that template)
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
    #  (b) if it DOES appear, compare react_adj[map_no] vs. prod_adj[map_no]
    #      If they differ -> bond‐pattern changed -> reactive
    #  (c) also check if atomic number or formal charge changed (rare in a template, but could)
    #      We compare the two atoms directly. We need to find the reactant‐template Atom and product‐template
    #      Atom to compare
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


def index_uncontested(
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

        reactive_inds = reactive_template_atoms(rl.rxn)[0]
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
    uncontested = [(rl, match_set) for rl, match_set in uncontested if (rl.id, frozenset(match_set)) not in failed_combos]

    return uncontested


def apply_uncontested(
    parent: Mol,
    uncontested: list[tuple[ReactionRule, set[int]]],
    original_taken_tags: set[int],
) -> tuple[list[Mol], list[ReactionRule], set[tuple[int, frozenset[int]]]]:
    """
    Apply uncontested reactions in bulk.

    :param parent: RDKit molecule
    :param uncontested: List of uncontested reactions
    :param original_taken_tags: List of atom tags from original reactant
    :return: List of trtue products, a list of applied ReactionRules,  and a set of failed combinations
    """
    originaL_taken_tags = list(original_taken_tags)

    applied_reactions: list[ReactionRule] = []

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

        unmasked_parent = copy_mol(parent)  # keep original parent for later
        results = rl(parent, msk)  # apply reaction rule

        try:
            if len(results) == 0:
                raise ValueError(f"No products from uncontested reaction {rl.rid}")

            if len(results) > 1:
                raise ValueError(f"More than one product from uncontested reaction {rl.rid}")

            result = results[0]
            applied_reactions.append(rl)  # keep track of successfully applied reactions

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
