"""Module contains functionality for reconstructing a linear backbone from RetroMol's linear readout."""

import logging
from dataclasses import dataclass
from enum import Enum
from typing import Any

from rdkit.Chem import Mol

from retromol.chem.mol import encode_mol, smiles_to_mol, smarts_to_mol, mol_to_smiles
from retromol.chem.tagging import get_tags_mol
from retromol.chem.reaction import smarts_to_reaction
from retromol.model.readout import LinearReadout
from retromol.model.result import Result

logger = logging.getLogger(__name__)


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
rxn_fuse_aa_alpha_aa_alpha = smarts_to_reaction(r"[N:1][C:2]C(=O)[OH].[N:3][C:4]-,=[C:5](=[O:6])[OH:7]>>[N:1][C:2][N:3][C:4]-,=[C:5](=[O:6])[OH:7]")


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


def react(rxn, reactants: tuple) -> Mol:
    """
    Apply an RDKit reaction to a set of reactants and return the first product.

    :param rxn: The RDKit reaction to apply.
    :param reactants: The reactant mols.
    :return: The first product mol of the first matched reactant combination.
    :raises BackboneReconstructionError: If the reaction produced no products.
    """
    products = rxn.RunReactants(reactants)
    if not products:
        raise BackboneReconstructionError(f"Reaction {rxn} produced no products for given reactants.")
    return products[0][0]


def fuse(mol: Mol, ext_mol: Mol, prev_type: MotifType | None, curr_type: MotifType | None) -> Mol:
    """
    Fuse two motifs into a single mol, using the fusion chemistry appropriate for
    the pair of motif types.

    :param mol: The mol built up so far.
    :param ext_mol: The next building-block mol to fuse onto it.
    :param prev_type: The motif type of the tail of `mol`.
    :param curr_type: The motif type of `ext_mol`.
    :return: The fused mol.
    :raises BackboneReconstructionError: If there is no fusion rule for this pair of motif types.
    """
    match prev_type, curr_type:
        case MotifType.STARTER, MotifType.PK_SINGLE | MotifType.PK_DOUBLE:
            return react(rxn_fuse_starter_pk, (mol, ext_mol))
        case MotifType.STARTER, MotifType.AA_ALPHA:
            return react(rxn_fuse_starter_aa_alpha, (mol, ext_mol))
        case MotifType.PK_SINGLE | MotifType.PK_DOUBLE, MotifType.PK_SINGLE | MotifType.PK_DOUBLE:
            return react(rxn_fuse_pk_pk, (mol, ext_mol))
        case MotifType.AA_ALPHA, MotifType.PK_SINGLE | MotifType.PK_DOUBLE:
            return react(rxn_fuse_aa_alpha_pk, (mol, ext_mol))
        case MotifType.PK_SINGLE | MotifType.PK_DOUBLE, MotifType.AA_ALPHA:
            return react(rxn_fuse_pk_aa_alpha, (mol, ext_mol))
        case MotifType.AA_ALPHA, MotifType.AA_ALPHA:
            return react(rxn_fuse_aa_alpha_aa_alpha, (mol, ext_mol))
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

    if starter is not None:
        prev_type = MotifType.STARTER
        prod = react(rxn_starter, (smiles_to_mol(starter),))

    for curr_smi in building_blocks:
        curr_mol = smiles_to_mol(curr_smi)
        curr_type = determine_type(curr_mol)

        if prod is None:
            if curr_type == MotifType.PK_SINGLE:
                prod = react(rxn_pk_start_single, (curr_mol,))
            elif curr_type == MotifType.PK_DOUBLE:
                prod = react(rxn_pk_start_double, (curr_mol,))
            else:
                prod = curr_mol
        else:
            if curr_type == MotifType.PK_SINGLE:
                curr_mol = react(rxn_pk_single, (curr_mol,))
                prod = fuse(prod, curr_mol, prev_type, curr_type)
            elif curr_type == MotifType.PK_DOUBLE:
                curr_mol = react(rxn_pk_double, (curr_mol,))
                prod = fuse(prod, curr_mol, prev_type, curr_type)
            elif curr_type == MotifType.AA_ALPHA:
                prod = fuse(prod, curr_mol, prev_type, curr_type)
            else:
                prod = curr_mol

        prev_type = curr_type

    if prod is None:
        raise BackboneReconstructionError("No backbone product could be built for this path.")

    return prod


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
            not (all(eligible[1:]) or all(eligible[:1]) or all(eligible))
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
            # Flip orientation
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
