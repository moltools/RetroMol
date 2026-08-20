"""Module contains functionality for reconstructing a linear backbone from RetroMol's linear readout."""

import logging
from dataclasses import dataclass
from enum import Enum
from typing import Any

from rdkit.Chem import Mol, RWMol
from rdkit.Chem.rdchem import Atom, BondStereo, BondType
from rdkit.Chem.rdmolops import AssignStereochemistry, SetDoubleBondNeighborDirections

from retromol.chem.mol import encode_mol, smiles_to_mol, smarts_to_mol, mol_to_smiles, sanitize_mol
from retromol.chem.tagging import get_tags_mol
from retromol.chem.reaction import smarts_to_reaction
from retromol.chem.stereo import BondStereoRecord
from retromol.model.readout import LinearReadout
from retromol.model.result import Result

logger = logging.getLogger(__name__)

_STEREO_TYPES = (BondStereo.STEREOE, BondStereo.STEREOZ, BondStereo.STEREOCIS, BondStereo.STEREOTRANS)
_CIS_TYPES = (BondStereo.STEREOZ, BondStereo.STEREOCIS)


def _capture_double_bond_stereo(mol: Mol) -> dict[frozenset[int], BondStereoRecord]:
    """
    Record the geometry of every stereo-defined double bond in `mol`, keyed by the
    isotope tags of its two atoms.

    A relaxed variant of `retromol.chem.stereo.capture_double_bond_stereo`: that
    function requires the stereo reference substituent itself to be tagged, which
    holds for the main retrobiosynthetic parsing pipeline (it captures stereo from
    the original input molecule, where every atom is tagged) but not here --
    reconstruction operates on already-parsed building-block fragments, where a
    double bond's reference substituent is often the (untagged) leaving-group
    placeholder atom introduced earlier during parsing, not an atom traceable to
    the original input compound. `restore_double_bond_stereo`'s own resolution logic
    already handles an untagged (0) reference correctly (it falls back to whichever
    single untagged heavy-atom substituent remains) -- only the capture step here
    needed relaxing to actually produce such an entry.

    :param mol: Molecule to inspect. Every atom is expected to already carry a
        unique isotope tag, or 0 if untagged (e.g. a leaving-group placeholder).
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

        registry[frozenset((tag_a, tag_b))] = BondStereoRecord(
            ref_tags={tag_a: ref_tag_a, tag_b: ref_tag_b},
            cis=stereo in _CIS_TYPES,
        )

    return registry


def _resolve_stereo_reference(atom: Atom, double_bond_partner: Atom, original_ref_tag: int) -> int | None:
    """
    Find the atom index that should now act as the stereo reference substituent on `atom`.

    A relaxed variant of `retromol.chem.stereo._resolve_stereo_reference`, with one
    extra fast path: if `atom` has exactly one substituent besides the double bond
    itself, that substituent is unambiguously the reference regardless of its tag --
    there's no other candidate it could be. The shared version doesn't need this,
    since it only ever looks up an `original_ref_tag` that was itself required to be
    a real tag; here `original_ref_tag` can be 0 (untagged in the original), and a
    single fusion step can legitimately turn that untagged placeholder into a *newly
    tagged* substituent in one move (a leaving-group placeholder atom being replaced
    by a bond straight to the existing, already-tagged backbone) -- a transition the
    shared logic's two cases (same tag persists, or persists as a fresh untagged
    atom) doesn't cover, but which is exactly what fusing two units together does.

    :param atom: One of the two double bond atoms.
    :param double_bond_partner: The atom on the other end of the double bond.
    :param original_ref_tag: Isotope tag of the substituent that used to define this
        side's geometry, before the reaction ran (0 if it was itself untagged).
    :return: Atom index of the resolved reference substituent, or None if it can't be
        determined unambiguously (in which case the caller should leave stereo unset).
    """
    substituents = [n for n in atom.GetNeighbors() if n.GetIdx() != double_bond_partner.GetIdx()]
    if not substituents:
        return None

    if len(substituents) == 1:
        return substituents[0].GetIdx()

    for n in substituents:
        if n.GetIsotope() == original_ref_tag:
            return n.GetIdx()

    new_substituents = [n for n in substituents if n.GetAtomicNum() > 1 and n.GetIsotope() == 0]
    if len(new_substituents) == 1:
        return new_substituents[0].GetIdx()

    return None


def _restore_double_bond_stereo(mol: Mol, registry: dict[frozenset[int], BondStereoRecord]) -> Mol:
    """
    Restore double bond stereo on `mol` for bonds that lost it during a reaction but
    whose two atoms match a stereo-defined bond recorded in `registry`.

    A copy of `retromol.chem.stereo.restore_double_bond_stereo` that resolves
    references via this module's relaxed `_resolve_stereo_reference` instead of the
    shared one -- see that function's docstring for why.

    :param mol: Product molecule to restore stereo on, mutated in place.
    :param registry: Registry produced by `_capture_double_bond_stereo`.
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
        SetDoubleBondNeighborDirections(mol)
        AssignStereochemistry(mol, cleanIt=True, force=True)

    return mol


class BackboneReconstructionError(RuntimeError):
    """Raised when the hardcoded fusion chemistry below cannot combine two motifs into a backbone."""


BACKBONE_WARNING = (
    "The linear backbone could not be reconstructed with RetroMol's current "
    "hardcoded fusion chemistry. The primary sequence below is still derived "
    "directly from the source structure, so individual units can still be "
    "highlighted there, but no reconstructed backbone structure is shown."
)

UNORDERED_WARNING = (
    "RetroMol identified all of these building blocks, but couldn't connect them "
    "into a single biosynthetic order (e.g. a branched or cyclic assembly), so no "
    "primary sequence could be assembled. They're shown below as an unordered set. "
    "Click one to highlight its atoms in the structure."
)


pattern_pk_single = smarts_to_mol(r"OSC-CC(=O)[OH]")
pattern_pk_double = smarts_to_mol("OSC=CC(=O)[OH]")
pattern_aa_alpha = smarts_to_mol(r"NCC(=O)[OH]")
pattern_aa_beta = smarts_to_mol(r"NCCC(=O)[OH]")
eligible_patterns = [pattern_pk_single, pattern_pk_double, pattern_aa_alpha, pattern_aa_beta]

rxn_starter = smarts_to_reaction(r"[*:1][C:2](=[O:3])[OH:4]>>[*:1][C:2](=[O:3])[O:4]-[SnH]")
rxn_pk_start_single = smarts_to_reaction(r"OS[C:1]-[C:2][C:3](=[O:4])[OH:5]>>[*][C:1]-[C:2][C:3](=[O:4])[O:5]-[SnH]")
rxn_pk_start_double = smarts_to_reaction(r"OS[C:1]=[C:2][C:3](=[O:4])[OH:5]>>[*][C:1]=[C:2][C:3](=[O:4])[O:5]-[SnH]")
rxn_pk_single = smarts_to_reaction(r"OS[C:1]-[C:2][C:3](=[O:4])[OH:5]>>[PbH][C:1]-[C:2][C:3](=[O:4])[O:5]-[SnH]")
rxn_pk_double = smarts_to_reaction(r"OS[C:1]=[C:2][C:3](=[O:4])[OH:5]>>[PbH][C:1]=[C:2][C:3](=[O:4])[O:5]-[SnH]")

rxn_fuse_starter_pk = smarts_to_reaction(r"[*:1]C(=O)O[SnH].[PbH][C:2]~[C:3][C:4](=[O:5])[O:6][SnH]>>[*:1][C:2]~[C:3][C:4](=[O:5])[O:6][SnH]")
rxn_fuse_starter_aa_alpha = smarts_to_reaction(r"[*:1][C:2](=[O:3])O[SnH].[N:4][C:5]C(=O)[OH]>>[*:1][C:2](=[O:3])[N:4][C:5]C(=O)[OH]")
rxn_fuse_pk_pk = smarts_to_reaction(r"[*:1][C:2]~[C:3]C(=O)O[SnH].[PbH][C:4]~[C:5][C:6](=[O:7])[O:8][SnH]>>[*:1][C:2]~[C:3][C:4]~[C:5][C:6](=[O:7])[O:8][SnH]")
rxn_fuse_aa_alpha_pk = smarts_to_reaction(r"[N:1][C:2]C(=O)[OH].[PbH][C:3]~[C:4][C:5](=[O:6])[O:7][SnH:8]>>[N:1][C:2][C:3]~[C:4][C:5](=[O:6])[O:7][SnH:8]")
rxn_fuse_pk_aa_alpha = smarts_to_reaction(r"[*:1][C:2]~[C:3][C:4](=[O:5])O[SnH].[N:6][C:7]C(=O)[OH]>>[*:1][C:2]~[C:3][C:4](=[O:5])[N:6][C:7]C(=O)[OH]")
rxn_fuse_aa_alpha_aa_alpha = smarts_to_reaction(r"[N:1][C:2][C:8](=[O:9])[OH].[N:3][C:4]-,=[C:5](=[O:6])[OH:7]>>[N:1][C:2][C:8](=[O:9])[N:3][C:4]-,=[C:5](=[O:6])[OH:7]")


class MotifType(Enum):
    STARTER = "STARTER"
    PK_SINGLE = "PK_SINGLE"
    PK_DOUBLE = "PK_DOUBLE"
    AA_ALPHA = "AA_ALPHA"


def determine_type(mol: Mol) -> MotifType | None:
    """
    Classify a building-block mol into one of the fusion-chemistry motif types below.

    :param mol: The building-block mol to classify.
    :return: The matched motif type, or None if it doesn't match any known pattern.
    """
    if mol.HasSubstructMatch(pattern_aa_alpha):
        return MotifType.AA_ALPHA
    elif mol.HasSubstructMatch(pattern_pk_single):
        return MotifType.PK_SINGLE
    elif mol.HasSubstructMatch(pattern_pk_double):
        return MotifType.PK_DOUBLE
    else:
        return None


def react(
    rxn,
    reactants: tuple,
    stereo_registry: dict[frozenset[int], BondStereoRecord] | None = None,
) -> Mol:
    """
    Apply an RDKit reaction to a set of reactants and return the first product.

    :param rxn: The RDKit reaction to apply.
    :param reactants: The reactant mols.
    :param stereo_registry: Registry of double bond stereo recorded from the original
        input molecule (see `capture_double_bond_stereo`), used to restore E/Z stereo
        on the product that this reaction's SMIRKS -- none of which encode `/`/`\\`
        stereo bond markers -- would otherwise silently discard. Updated in place with
        any stereo newly present on the product, so later reactions in the same chain
        can restore it too. Skipped when None.
    :return: The first product mol of the first matched reactant combination.
    :raises BackboneReconstructionError: If the reaction produced no products.
    """
    products = rxn.RunReactants(reactants)
    if not products:
        raise BackboneReconstructionError(f"Reaction {rxn} produced no products for given reactants.")
    product = products[0][0]

    if stereo_registry is not None:
        if stereo_registry:
            _restore_double_bond_stereo(product, stereo_registry)
        for key, record in _capture_double_bond_stereo(product).items():
            stereo_registry.setdefault(key, record)

    return product


def fuse(
    mol: Mol,
    ext_mol: Mol,
    prev_type: MotifType | None,
    curr_type: MotifType | None,
    stereo_registry: dict[frozenset[int], BondStereoRecord] | None = None,
) -> Mol:
    """
    Fuse two motifs into a single mol, using the fusion chemistry appropriate for
    the pair of motif types.

    :param mol: The mol built up so far.
    :param ext_mol: The next building-block mol to fuse onto it.
    :param prev_type: The motif type of the tail of `mol`.
    :param curr_type: The motif type of `ext_mol`.
    :param stereo_registry: Forwarded to `react` -- see there.
    :return: The fused mol.
    :raises BackboneReconstructionError: If there is no fusion rule for this pair of motif types.
    """
    match prev_type, curr_type:
        case MotifType.STARTER, MotifType.PK_SINGLE | MotifType.PK_DOUBLE:
            return react(rxn_fuse_starter_pk, (mol, ext_mol), stereo_registry)
        case MotifType.STARTER, MotifType.AA_ALPHA:
            return react(rxn_fuse_starter_aa_alpha, (mol, ext_mol), stereo_registry)
        case MotifType.PK_SINGLE | MotifType.PK_DOUBLE, MotifType.PK_SINGLE | MotifType.PK_DOUBLE:
            return react(rxn_fuse_pk_pk, (mol, ext_mol), stereo_registry)
        case MotifType.AA_ALPHA, MotifType.PK_SINGLE | MotifType.PK_DOUBLE:
            return react(rxn_fuse_aa_alpha_pk, (mol, ext_mol), stereo_registry)
        case MotifType.PK_SINGLE | MotifType.PK_DOUBLE, MotifType.AA_ALPHA:
            return react(rxn_fuse_pk_aa_alpha, (mol, ext_mol), stereo_registry)
        case MotifType.AA_ALPHA, MotifType.AA_ALPHA:
            return react(rxn_fuse_aa_alpha_aa_alpha, (mol, ext_mol), stereo_registry)
        case _:
            raise BackboneReconstructionError(f"No fusion rule for motif transition {prev_type} -> {curr_type}.")


@dataclass(frozen=True)
class Reconstruction:
    """
    Reconstruct a linear backbone from RetroMol's linear readout.

    :param tagged_input_smiles: Input SMILES string with tagged atoms.
    :param tagged_backbone_smiles: Reconstructed backbone SMILES string with tagged atoms,
        corresponding to tagged_input_smiles, or None if the backbone could not be
        reconstructed (see `backbone_warning` in that case).
    :param primary_sequence: Per-module name and tags for items present in the primary sequence.
    :param backbone_warning: Set when `tagged_backbone_smiles` is None, explaining why the
        backbone is missing even though a primary sequence was found.
    :param ordered: False when `primary_sequence` isn't a genuine biosynthetic order --
        RetroMol identified every item individually, but couldn't connect them into a
        single path (e.g. a branched or cyclic assembly), so they're just the set of
        building blocks found, in no particular order.
    """

    tagged_input_smiles: str
    tagged_backbone_smiles: str | None
    primary_sequence: list[tuple[str, set[int]]]
    backbone_warning: str | None = None
    ordered: bool = True

    def to_dict(self) -> dict[str, Any]:
        """
        Convert the reconstruction to a dictionary.

        :return: A dictionary representation of the reconstruction.
        """
        return {
            "tagged_input_smiles": self.tagged_input_smiles,
            "tagged_backbone_smiles": self.tagged_backbone_smiles,
            "primary_sequence": [(x, list(y)) for x, y in self.primary_sequence],
            "backbone_warning": self.backbone_warning,
            "ordered": self.ordered,
        }


def _reconstruct_backbone(starter: str | None, building_blocks: list[str]) -> Mol:
    """
    Build the linear backbone product mol for one ordered list of building blocks.

    :param starter: SMILES of the non-eligible starter unit, or None if every unit
        in the path is itself an eligible building block.
    :param building_blocks: Ordered building-block SMILES to fuse onto the starter.
    :return: The fused backbone mol.
    :raises BackboneReconstructionError: If any step of the fusion chemistry fails.
    """
    prev_type: MotifType | None = None
    prod: Mol | None = None

    # Recorded once per unit, from its tagged mol as originally identified in the
    # compound -- *before* any fusion reaction runs on it -- then threaded through
    # every reaction below so `react` can restore E/Z stereo the fusion SMIRKS
    # themselves don't encode (see `react`'s docstring).
    stereo_registry: dict[frozenset[int], BondStereoRecord] = {}

    if starter is not None:
        starter_mol = smiles_to_mol(starter)
        for key, record in _capture_double_bond_stereo(starter_mol).items():
            stereo_registry.setdefault(key, record)
        prev_type = MotifType.STARTER
        prod = react(rxn_starter, (starter_mol,), stereo_registry)

    for curr_smi in building_blocks:
        curr_mol = smiles_to_mol(curr_smi)
        for key, record in _capture_double_bond_stereo(curr_mol).items():
            stereo_registry.setdefault(key, record)
        curr_type = determine_type(curr_mol)

        if prod is None:
            if curr_type == MotifType.PK_SINGLE:
                prod = react(rxn_pk_start_single, (curr_mol,), stereo_registry)
            elif curr_type == MotifType.PK_DOUBLE:
                prod = react(rxn_pk_start_double, (curr_mol,), stereo_registry)
            else:
                prod = curr_mol
        else:
            if curr_type == MotifType.PK_SINGLE:
                curr_mol = react(rxn_pk_single, (curr_mol,), stereo_registry)
                prod = fuse(prod, curr_mol, prev_type, curr_type, stereo_registry)
            elif curr_type == MotifType.PK_DOUBLE:
                curr_mol = react(rxn_pk_double, (curr_mol,), stereo_registry)
                prod = fuse(prod, curr_mol, prev_type, curr_type, stereo_registry)
            elif curr_type == MotifType.AA_ALPHA:
                prod = fuse(prod, curr_mol, prev_type, curr_type, stereo_registry)
            else:
                # `curr_type` is None: this building block didn't match any known
                # fusion motif type (e.g. a beta-amino acid, which `eligible_patterns`
                # allows through but `determine_type` doesn't classify). There's no
                # fusion chemistry for it -- fail loudly rather than silently
                # discarding everything fused so far and replacing it with just this
                # unit, which would report a confidently wrong backbone as if it had
                # succeeded.
                raise BackboneReconstructionError(
                    "Could not classify a mid-chain building block for fusion (no matching motif type)."
                )

        prev_type = curr_type

    if prod is None:
        raise BackboneReconstructionError("No backbone product could be built for this path.")

    return _cap_dangling_placeholders(prod)


def _cap_dangling_placeholders(mol: Mol) -> Mol:
    """
    Replace any leftover leaving-group/activation placeholder atoms with a plain
    hydrogen, so a successfully-built backbone doesn't visibly contain fake atoms
    that were only ever meant as internal bonding handles during fusion:

    - `[SnH]` (always as `-O-[SnH]`): marks a unit's tail as "still activated,
      ready for the next condensation" -- if nothing ever fused onto it (the chain
      simply ends there), the natural resting state is a plain carboxylic acid, so
      it's replaced with `-OH`.
    - `[PbH]`: marks a unit's alpha carbon as "expects an incoming bond" --
      shouldn't normally survive to a *successfully* fused backbone, but is capped
      defensively the same way if it does.
    - An unfilled wildcard `*`: left at the very front of a chain that has no
      explicit non-eligible starter unit (see `rxn_pk_start_single`/`_double`).

    :param mol: The fused backbone mol, possibly still carrying placeholder atoms.
    :return: The same structure with placeholders replaced by an implicit H.
    :raises BackboneReconstructionError: If capping leaves the structure unsanitizable.
    """
    rw = RWMol(mol)
    to_remove: list[int] = []
    neighbors_to_reset: list[int] = []

    for atom in rw.GetAtoms():
        if atom.GetSymbol() in ("Sn", "Pb") or atom.GetAtomicNum() == 0:
            if atom.GetDegree() != 1:
                continue
            to_remove.append(atom.GetIdx())
            neighbors_to_reset.append(atom.GetNeighbors()[0].GetIdx())

    if not to_remove:
        return mol

    for idx in neighbors_to_reset:
        neighbor = rw.GetAtomWithIdx(idx)
        # The placeholder's neighbor was sanitized alongside it with its H count
        # frozen to accommodate the placeholder's own unusual valence -- reset it
        # so removing the placeholder lets RDKit recompute a normal H count instead
        # of leaving an under-valent radical.
        neighbor.SetNoImplicit(False)
        neighbor.SetNumExplicitHs(0)
        neighbor.SetNumRadicalElectrons(0)

    for idx in sorted(to_remove, reverse=True):
        rw.RemoveAtom(idx)

    capped = rw.GetMol()
    if not sanitize_mol(capped, fix_hydrogens=True):
        raise BackboneReconstructionError("Failed to cap dangling placeholder atoms in the assembled backbone.")

    return capped


def reconstruct_linear_readout(result: Result) -> list[Reconstruction]:
    """
    Reconstruct primary sequences (and, where possible, a linear backbone) from
    RetroMol's linear readout.

    A candidate path can yield a primary sequence without a reconstructed backbone:
    the backbone is rebuilt with a small hardcoded set of fusion reactions that
    doesn't cover every motif combination, and any failure there is reported via
    `Reconstruction.backbone_warning` rather than dropping the whole candidate.

    If no path could be threaded into a sequence at all (e.g. the molecule's
    building blocks don't form a single linear chain), but building blocks were
    still individually identified, a single unordered `Reconstruction`
    (`ordered=False`) listing them is returned instead of an empty list -- see
    `UNORDERED_WARNING`.

    :param result: The result object returned by RetroMol, containing the reaction graph and the original molecule.
    :return: The reconstructed candidates, one per eligible path. Empty if no path
        in the readout consists of (mostly) eligible primary-sequence building blocks.
    """
    root_enc = encode_mol(result.submission.mol)
    readout = LinearReadout.from_reaction_graph(root_enc, reaction_graph=result.reaction_graph, identified_only=True)
    tagged_input_smiles = mol_to_smiles(result.submission.mol, include_tags=True)

    reconstructions: list[Reconstruction] = []

    paths = sorted(readout.paths, key=lambda x: len(x), reverse=True)
    for path in paths:
        building_blocks: list[str] = []
        primary_sequence: list[tuple[str, set[int]]] = []
        eligible = [False for _ in range(len(path))]

        for item_idx, item in enumerate(path):
            item_name: str = item.identity.matched_rule.name if item.identified else "X"
            item_tags: set[int] = get_tags_mol(item.mol)
            item_smiles = mol_to_smiles(item.mol, include_tags=True)

            if any(item.mol.HasSubstructMatch(pattern) for pattern in eligible_patterns):
                eligible[item_idx] = True

            primary_sequence.append((item_name, item_tags))
            building_blocks.append(item_smiles)

        # An eligible path has all items as eligible True, but can have first or last item as False (this is the starting unit)
        # Path also has to be at least 2 units long
        if (
            not (all(eligible[1:]) or all(eligible[:-1]) or all(eligible))
            or len(eligible) <= 1
        ):
            continue

        # Re-orient path to have non-eligible item at start of sequence
        if all(eligible):
            # Nothing to orient, all eligible building blocks
            starter = None
        elif all(eligible[1:]):
            # Already good orientation, but remove the starter
            starter = building_blocks[0]
            building_blocks = building_blocks[1:]
        else:
            # all(eligible[:-1]) -- non-eligible starter is at the *end*, flip so it's first
            building_blocks = list(reversed(building_blocks))
            primary_sequence = list(reversed(primary_sequence))

            # Remove starter
            starter = building_blocks[0]
            building_blocks = building_blocks[1:]

        if not building_blocks:
            continue

        tagged_backbone_smiles: str | None = None
        backbone_warning: str | None = None
        try:
            backbone_mol = _reconstruct_backbone(starter, building_blocks)
            tagged_backbone_smiles = mol_to_smiles(backbone_mol, include_tags=True)
        except Exception:
            logger.warning(
                "reconstruct_linear_readout: backbone reconstruction failed for a candidate path",
                exc_info=True,
            )
            backbone_warning = BACKBONE_WARNING

        reconstructions.append(
            Reconstruction(
                tagged_input_smiles=tagged_input_smiles,
                tagged_backbone_smiles=tagged_backbone_smiles,
                primary_sequence=primary_sequence,
                backbone_warning=backbone_warning,
            )
        )

    # No path threaded its identified building blocks into a usable primary
    # sequence (e.g. a branched or cyclic assembly RetroMol's path-finding
    # doesn't linearize), but individual building blocks may still have been
    # identified across the molecule. Surface those as an explicitly unordered
    # set rather than reporting nothing, so a fully (or mostly) covered compound
    # doesn't look unparsed just because no single order could be determined.
    if not reconstructions:
        monomer_nodes = readout.assembly_graph.monomer_nodes()
        if monomer_nodes:
            unordered_sequence = [
                (node.identity.matched_rule.name if node.identified else "X", get_tags_mol(node.mol))
                for node in monomer_nodes
            ]
            reconstructions.append(
                Reconstruction(
                    tagged_input_smiles=tagged_input_smiles,
                    tagged_backbone_smiles=None,
                    primary_sequence=unordered_sequence,
                    backbone_warning=UNORDERED_WARNING,
                    ordered=False,
                )
            )

    return reconstructions
