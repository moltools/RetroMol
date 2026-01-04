"""Module for RDKit reaction utilities."""

from rdkit.Chem.rdChemReactions import ChemicalReaction, ReactionFromSmarts


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
