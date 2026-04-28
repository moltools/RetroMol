"""Module contains functionality for reconstructing a linear backbone from RetroMol's linear readout."""

import re
from dataclasses import dataclass
from enum import Enum
from typing import Any

from retromol.chem.mol import encode_mol, smiles_to_mol, smarts_to_mol, mol_to_smiles
from retromol.chem.tagging import get_tags_mol
from retromol.chem.reaction import smarts_to_reaction
from retromol.model.readout import LinearReadout
from retromol.model.result import Result


pattern_pk_single = smarts_to_mol(r"OSC-CC(=O)[OH]")
pattern_pk_double = smarts_to_mol("OSC=CC(=O)[OH]")
pattern_aa_alpha = smarts_to_mol(r"NCC(=O)[OH]")
pattern_aa_beta = smarts_to_mol(r"NCCC(=O)[OH]")
eligible_patterns = [pattern_pk_single, pattern_pk_double, pattern_aa_alpha, pattern_aa_beta]

# rxn_starter = smarts_to_reaction(r"[*:1][C](=O)[OH]>>[*:1][SnH]")
rxn_pk_start_single = smarts_to_reaction(r"OS[C:1]-[C:2][C:3](=[O:4])[OH:5]>>[*][C:1]-[C:2][C:3](=[O:4])[O:5]-[SnH]")
rxn_pk_start_double = smarts_to_reaction(r"OS[C:1]=[C:2][C:3](=[O:4])[OH:5]>>[*][C:1]=[C:2][C:3](=[O:4])[O:5]-[SnH]")
rxn_pk_single = smarts_to_reaction(r"OS[C:1]-[C:2][C:3](=[O:4])[OH:5]>>[PbH][C:1]-[C:2][C:3](=[O:4])[O:5]-[SnH]")
rxn_pk_double = smarts_to_reaction(r"OS[C:1]=[C:2][C:3](=[O:4])[OH:5]>>[PbH][C:1]=[C:2][C:3](=[O:4])[O:5]-[SnH]")

rxn_fuse_pk_pk = smarts_to_reaction(r"[*:1][C:2]~[C:3]C(=O)O[SnH].[PbH][C:4]~[C:5][C:6](=[O:7])[O:8][SnH]>>[*:1][C:2]~[C:3][C:4]~[C:5][C:6](=[O:7])[O:8][SnH]")
rxn_fuse_aa_alpha_pk = smarts_to_reaction(r"[N:1][C:2]C(=O)[OH].[PbH][C:3]~[C:4][C:5](=[O:6])[O:7][SnH:8]>>[N:1][C:2][C:3]~[C:4][C:5](=[O:6])[O:7][SnH:8]")
rxn_fuse_pk_aa_alpha = smarts_to_reaction(r"[*:1][C:2]~[C:3]C(=O)O[SnH].[N:4][C:5]C(=O)[OH]>>[*:1][C:2]~[C:3][N:4][C:5]C(=O)[OH]")
rxn_fuse_aa_alpha_aa_alpha = smarts_to_reaction(r"[N:1][C:2]C(=O)[OH].[N:3][C:4]-,=[C:5](=[O:6])[OH:7]>>[N:1][C:2][N:3][C:4]-,=[C:5](=[O:6])[OH:7]")


class MotifType(Enum):
    PK_SINGLE = "PK_SINGLE"
    PK_DOUBLE = "PK_DOUBLE"
    AA_ALPHA = "AA_ALPHA"
    AA_BETA = "AA_BETA"


def determine_type(mol) -> MotifType | None:
    print("determining...", mol_to_smiles(mol))
    if mol.HasSubstructMatch(pattern_aa_alpha):
        print("AA alpha")
        return MotifType.AA_ALPHA
    elif mol.HasSubstructMatch(pattern_aa_beta):
        print("AA beta")
        return MotifType.AA_BETA
    elif mol.HasSubstructMatch(pattern_pk_single):
        print("PK")
        return MotifType.PK_SINGLE
    elif mol.HasSubstructMatch(pattern_pk_double):
        print("PK double")
        return MotifType.PK_DOUBLE
    else:
        return None


def fuse(mol, ext_mol, prev_type, curr_type):
    match prev_type, curr_type:
        case MotifType.PK_SINGLE | MotifType.PK_DOUBLE, MotifType.PK_SINGLE | MotifType.PK_DOUBLE:
            print("FUSE PK-PK", prev_type, curr_type)
            print(mol_to_smiles(mol), mol_to_smiles(ext_mol))
            return react(rxn_fuse_pk_pk, (mol, ext_mol,))
        case MotifType.AA_ALPHA, MotifType.PK_SINGLE | MotifType.PK_DOUBLE:
            print("FUSE AA alpha-PK", prev_type, curr_type)
            print(mol_to_smiles(mol), mol_to_smiles(ext_mol))
            return react(rxn_fuse_aa_alpha_pk, (mol, ext_mol,))
        case MotifType.PK_SINGLE | MotifType.PK_DOUBLE, MotifType.AA_ALPHA:
            print("FUSE PK-AA alpha", prev_type, curr_type)
            print(mol_to_smiles(mol), mol_to_smiles(ext_mol))
            return react(rxn_fuse_pk_aa_alpha, (mol, ext_mol,))
        case MotifType.AA_ALPHA, MotifType.AA_ALPHA:
            print("FUSE AA alpha-AA alpha", prev_type, curr_type)
            print(mol_to_smiles(mol), mol_to_smiles(ext_mol))
            return react(rxn_fuse_aa_alpha_aa_alpha, (mol, ext_mol,))
        case _:
            print("DIFF")
            print(prev_type, curr_type)
            return mol


def react(rxn, reactants):
    products = rxn.RunReactants(reactants)
    print(products)
    product = products[0][0]
    return product


@dataclass(frozen=True)
class Reconstruction:
    """
    Reconstruct a linear backbone from RetroMol's linear readout.

    :param tagged_input_smiles: Input SMILES string with tagged atoms.
    :param tagged_backbone_smiles: Reconstructed backbone SMILES string with tagged atoms, corresponding to tagged_input_smiles.
    :param primary_sequence: Per-module name and tags for items present in backbone SMILES.
    """

    tagged_input_smiles: str
    tagged_backbone_smiles: str
    primary_sequence: list[tuple[str, set[int]]]

    def to_dict(self) -> dict[str, Any]:
        """
        Convert the reconstruction to a dictionary.

        :return: A dictionary representation of the reconstruction.
        """
        return {
            "tagged_input_smiles": self.tagged_input_smiles,
            "tagged_backbone_smiles": self.tagged_backbone_smiles,
            "primary_sequence": [(x, list(y)) for x, y in self.primary_sequence],
        }


def reconstruct_linear_readout(result: Result) -> list[Reconstruction]:
    """
    Reconstruct a linear backbone from RetroMol's linear readout.

    :param result: The result object returned by RetroMol, containing the reaction graph and the original molecule.
    :return: The reconstructed backbone.
    :raises ValueError: If no successfull reconstructions are found.
    """
    root_enc = encode_mol(result.submission.mol)
    readout = LinearReadout.from_reaction_graph(root_enc, reaction_graph=result.reaction_graph, identified_only=True)

    reconstructions: list[Reconstruction] = []

    # Get building blocks
    paths = readout.paths
    paths.sort(key=lambda x: len(x), reverse=True)
    for path in readout.paths:

        building_blocks = []
        primary_sequence: list[tuple[str, set[int]]] = []

        eligible = [False for _ in range(len(path))]
        for item_idx, item in enumerate(path):
            item_name: str = item.identity.matched_rule.name if item.identified else "X"
            item_tags: set[int] = get_tags_mol(item.mol)
            item_smiles = mol_to_smiles(item.mol, include_tags=True)

            # For every item in path check if it is an eligible primary sequence building block
            if any([item.mol.HasSubstructMatch(pattern) for pattern in eligible_patterns]):
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

        print(eligible)
        # Re-orient path to have non-eligible item at start of sequence
        if all(eligible):
            # Nothing to orient, all aligible building blocks
            pass
        elif all(eligible[1:]):
            # Already good orientation, but remove the starter
            eligible = eligible[1:]
            building_blocks = building_blocks[1:]
            primary_sequence = primary_sequence[1:]
        else:
            # Flip orientation
            eligible.reverse()
            building_blocks.reverse()
            primary_sequence.reverse()

            # Remove starter
            eligible = eligible[1:]
            building_blocks = building_blocks[1:]
            primary_sequence = primary_sequence[1:]

        # Reconstruct linear backbone, we can have either PK-PK, PK-AA, AA-AA or AA-AA
        if not len(building_blocks):
            continue
        prev_type = None
        prod = None
        while building_blocks:
            curr_smi = building_blocks.pop(0)
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
            continue

        # If product contains Sn atom, replace it with wildcard
        print(mol_to_smiles(prod))

        product_smi = mol_to_smiles(prod, include_tags=True)
        input_smi = mol_to_smiles(result.submission.mol, include_tags=True)

        reconstruction = Reconstruction(
            tagged_input_smiles=input_smi,
            tagged_backbone_smiles=product_smi,
            primary_sequence=primary_sequence,
        )
        reconstructions.append(reconstruction)

    if len(reconstructions) == 0:
        raise ValueError("No successfull reconstructions!")

    return reconstructions
