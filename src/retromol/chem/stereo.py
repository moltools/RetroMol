"""Module for restoring double bond (E/Z) stereochemistry that reaction rules discard.

Reaction SMARTS that replace one of a stereo double bond's substituents (e.g. cleaving
a chain carbon off an alkene to attach a thiol leaving-group placeholder) leave RDKit
with no information to derive the product bond's stereo from -- the substituent that
defined the geometry is gone, and reaction product templates in `rxn.yml` don't encode
`/`/`\\` themselves. Since the atoms of the double bond usually keep their original
isotope tags across a reaction, we can look up what stereo *that exact bond* (identified
by its two atom tags) had in the original input molecule and, when the replacement
substituent occupies the same graph position 1:1, carry the geometry forward.

This is a no-op for bonds that never existed as a stereo-defined double bond in the
original input (e.g. a bond created by "shifted" reversed-polyketide-synthesis rules,
where the product's double bond sits at a different atom pair than any double bond in
the input) -- there is nothing to look up for those, and none is invented.
"""

from dataclasses import dataclass

from rdkit.Chem.rdchem import Atom, Bond, BondStereo, BondType, Mol
from rdkit.Chem.rdmolops import AssignStereochemistry, SetDoubleBondNeighborDirections

_STEREO_TYPES = (BondStereo.STEREOE, BondStereo.STEREOZ, BondStereo.STEREOCIS, BondStereo.STEREOTRANS)
_CIS_TYPES = (BondStereo.STEREOZ, BondStereo.STEREOCIS)


@dataclass(frozen=True)
class BondStereoRecord:
    """
    Recorded geometry of one stereo-defined double bond, keyed elsewhere by the isotope
    tags of its two atoms.

    :cvar ref_tags: Maps each of the bond's two atom tags to the isotope tag of the
        neighbor atom that was used as its stereo reference substituent.
    :cvar cis: Whether the two reference substituents are cis (True) or trans (False)
        to each other.
    """

    ref_tags: dict[int, int]
    cis: bool


def capture_double_bond_stereo(mol: Mol) -> dict[frozenset[int], BondStereoRecord]:
    """
    Record the geometry of every stereo-defined double bond in `mol`, keyed by the
    isotope tags of its two atoms.

    :param mol: Molecule to inspect. Every atom is expected to already carry a unique,
        nonzero isotope tag (see `tag_mol`).
    :return: Registry of double bond stereo records, keyed by the frozenset of the
        bond's two atom tags.
    """
    registry: dict[frozenset[int], BondStereoRecord] = {}

    for bond in mol.GetBonds():
        if bond.GetBondType() != BondType.DOUBLE:
            continue

        stereo = bond.GetStereo()
        if stereo not in _STEREO_TYPES:
            continue

        stereo_atom_idxs = bond.GetStereoAtoms()
        if len(stereo_atom_idxs) != 2:
            continue

        a, b = bond.GetBeginAtom(), bond.GetEndAtom()
        tag_a, tag_b = a.GetIsotope(), b.GetIsotope()
        if tag_a == 0 or tag_b == 0 or tag_a == tag_b:
            continue

        ref_0 = mol.GetAtomWithIdx(stereo_atom_idxs[0])
        ref_1 = mol.GetAtomWithIdx(stereo_atom_idxs[1])
        a_neighbor_idxs = {n.GetIdx() for n in a.GetNeighbors()}
        b_neighbor_idxs = {n.GetIdx() for n in b.GetNeighbors()}

        if ref_0.GetIdx() in a_neighbor_idxs and ref_1.GetIdx() in b_neighbor_idxs:
            ref_tag_a, ref_tag_b = ref_0.GetIsotope(), ref_1.GetIsotope()
        elif ref_0.GetIdx() in b_neighbor_idxs and ref_1.GetIdx() in a_neighbor_idxs:
            ref_tag_a, ref_tag_b = ref_1.GetIsotope(), ref_0.GetIsotope()
        else:
            continue

        if ref_tag_a == 0 or ref_tag_b == 0:
            continue

        registry[frozenset((tag_a, tag_b))] = BondStereoRecord(
            ref_tags={tag_a: ref_tag_a, tag_b: ref_tag_b},
            cis=stereo in _CIS_TYPES,
        )

    return registry


def _resolve_stereo_reference(atom: Atom, double_bond_partner: Atom, original_ref_tag: int) -> int | None:
    """
    Find the atom index that should now act as the stereo reference substituent on `atom`.

    :param atom: One of the two double bond atoms.
    :param double_bond_partner: The atom on the other end of the double bond.
    :param original_ref_tag: Isotope tag of the substituent that used to define this
        side's geometry, before the reaction ran.
    :return: Atom index of the resolved reference substituent, or None if it can't be
        determined unambiguously (in which case the caller should leave stereo unset).
    """
    substituents = [n for n in atom.GetNeighbors() if n.GetIdx() != double_bond_partner.GetIdx()]
    if not substituents:
        return None

    # The original reference substituent is still directly attached: reuse it as-is.
    for n in substituents:
        if n.GetIsotope() == original_ref_tag:
            return n.GetIdx()

    # It was replaced (e.g. by a new leaving-group atom introduced by the reaction). Any
    # other substituent still carrying its own original tag (this matters for
    # trisubstituted alkenes, which have such a second, untouched substituent alongside
    # the reference one) is not a candidate -- it was never the reference. Only a
    # substituent freshly introduced by this same reaction (no tag yet) can be the
    # replacement; if there's exactly one, the reaction swapped it in 1:1 at the same
    # graph position, so the geometric relationship still holds through it.
    new_substituents = [n for n in substituents if n.GetAtomicNum() > 1 and n.GetIsotope() == 0]
    if len(new_substituents) == 1:
        return new_substituents[0].GetIdx()

    return None


def restore_double_bond_stereo(mol: Mol, registry: dict[frozenset[int], BondStereoRecord]) -> Mol:
    """
    Restore double bond stereo on `mol` for bonds that lost it during a reaction but
    whose two atoms match a stereo-defined bond recorded in `registry`.

    :param mol: Product molecule to restore stereo on, mutated in place.
    :param registry: Registry produced by `capture_double_bond_stereo` on the original
        input molecule.
    :return: The same molecule, for convenient chaining.
    .. note:: This function mutates the input molecule in place.
    """
    if not registry:
        return mol

    restored_any = False

    for bond in mol.GetBonds():
        if bond.GetBondType() != BondType.DOUBLE or bond.GetStereo() != BondStereo.STEREONONE:
            continue

        a, b = bond.GetBeginAtom(), bond.GetEndAtom()
        tag_a, tag_b = a.GetIsotope(), b.GetIsotope()
        if tag_a == 0 or tag_b == 0:
            continue

        record = registry.get(frozenset((tag_a, tag_b)))
        if record is None:
            continue

        ref_a_idx = _resolve_stereo_reference(a, b, record.ref_tags[tag_a])
        ref_b_idx = _resolve_stereo_reference(b, a, record.ref_tags[tag_b])
        if ref_a_idx is None or ref_b_idx is None:
            continue

        bond.SetStereoAtoms(ref_a_idx, ref_b_idx)
        bond.SetStereo(BondStereo.STEREOCIS if record.cis else BondStereo.STEREOTRANS)
        restored_any = True

    if restored_any:
        # `/`/`\` SMILES output is driven by BondDir on the neighboring single bonds,
        # not directly by the double bond's Stereo/StereoAtoms -- derive it from what
        # we just set, then re-perceive so the assignment is fully consistent.
        SetDoubleBondNeighborDirections(mol)
        AssignStereochemistry(mol, cleanIt=True, force=True)

    return mol
