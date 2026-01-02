"""Module for masking atoms in a molecule based on their tags."""

from rdkit.Chem.rdchem import Mol

from retromol.chem.tagging import all_atoms_have_unique_tags


def mask_atoms(mol: Mol, mask_tags: set[int]) -> dict[int, int]:
    """
    Set atom numbers not in mask to 0 based on atom tags.

    :param mol: the molecule to mask
    :param mask_tags: mask of atom tags
    :return: mapping of atom tags to pre-mask atomic numbers
    """
    # Ensure all atoms have unique tags
    if not all_atoms_have_unique_tags(mol):
        raise ValueError("all atoms in the molecule must have unique tags before masking")
    
    # Create mapping of tags to atomic numbers
    tag_to_atomic_num = {a.GetIsotope(): a.GetAtomicNum() for a in mol.GetAtoms()}

    # Now we can apply the mask, set atomic numbers to 0 for atoms that do not have a tag in mask_tags
    for atom in mol.GetAtoms():
        if atom.GetIsotope() not in mask_tags:
            atom.SetAtomicNum(0)

    return tag_to_atomic_num


def unmask_atoms(mol: Mol, tag_to_atomic_num: dict[int, int]) -> None:
    """
    Restore atomic numbers based on the provided tag to atomic number mapping.

    :param mol: the molecule to unmask
    :param tag_to_atomic_num: mapping of atom tags to original atomic numbers
    :return: None
    .. note:: this function modifies the input molecule in place
    """
    for atom in mol.GetAtoms():
        if atom.GetIsotope() != 0:  # newly added atoms by rxn have isotope 0
            original_atom_num = tag_to_atomic_num.get(atom.GetIsotope(), None)
            if original_atom_num is None:
                raise ValueError(f"tag {atom.GetIsotope()} not found in tag_to_atomic_num mapping during unmasking")

            atom.SetAtomicNum(original_atom_num)

    # Validate that no atomic numbers are zero after unmasking
    curr_atomic_nums = {atom.GetAtomicNum() for atom in mol.GetAtoms()}
    if 0 in curr_atomic_nums:
        raise ValueError("unmasking failed, some atomic numbers are still zero")
    

def mapped_tags_changed(reactant: Mol, product: Mol) -> set[int]:
    """
    Determine which mapped atom tags have changed between reactant and product.

    :param reactant: reactant molecule
    :param product: product molecule
    :return: set of changed atom tags
    .. note:: a tag is considered changed if either the atom with that tag changes 
        atomic number, or its neighbor signature (by tagged IDs / atom types + bond order)
        changes
    """
    def _tag_to_idx(m: Mol) -> dict[int, int]:
        # Atom tags live in the isotope numbers; ignore zeros
        d: dict[int, int] = {}
        for a in m.GetAtoms():
            t = a.GetIsotope()
            if t:
                d[t] = a.GetIdx()
        return d

    def _neighbor_sig(m: Mol, ai: int) -> list[tuple[int, float]]:
        # Neighbor signature by (neighbor tag or neighbor atomicnum if untagged, bond order)
        out: list[tuple[int, float]] = []
        a = m.GetAtomWithIdx(ai)
        for b in a.GetBonds():
            nb = b.GetOtherAtomIdx(ai)
            na = m.GetAtomWithIdx(nb)
            ntag = na.GetIsotope()
            key = ntag if ntag else -na.GetAtomicNum()
            out.append((key, float(b.GetBondTypeAsDouble())))
        out.sort()
        return out

    changed: set[int] = set()
    rmap = _tag_to_idx(reactant)
    pmap = _tag_to_idx(product)
    for t in set(rmap).intersection(pmap):
        ra = reactant.GetAtomWithIdx(rmap[t])
        pa = product.GetAtomWithIdx(pmap[t])
        if ra.GetAtomicNum() != pa.GetAtomicNum():
            changed.add(t)
            continue
        if _neighbor_sig(reactant, rmap[t]) != _neighbor_sig(product, pmap[t]):
            changed.add(t)

    return changed


def is_masked_preserved(reactant: Mol, products: list[Mol], allowed: set[int]) -> bool:
    """
    Check that only tags in 'allowed' are changed across all products.

    :param reactant: reactant molecule
    :param products: list of product molecules
    :param allowed: set of allowed atom tags that can change
    :return: True if only allowed tags are changed, False otherwise
    """
    if not allowed:
        return True  # there is no mask, everything is allowed
    
    changed: set[int] = set()
    for pr in products:
        changed |= mapped_tags_changed(reactant, pr)

    return changed.issubset(allowed)
