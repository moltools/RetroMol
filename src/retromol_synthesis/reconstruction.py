"""Module contains functionality for reconstructing a linear backbone from RetroMol's linear readout."""

from dataclasses import dataclass
from typing import Any

from retromol.chem.mol import encode_mol, smiles_to_mol, smarts_to_mol, mol_to_smiles
from retromol.chem.tagging import get_tags_mol
from retromol.chem.reaction import smarts_to_reaction
from retromol.model.readout import LinearReadout
from retromol.model.result import Result


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

        print(len(path))

        building_blocks = []
        primary_sequence: list[tuple[str, set[int]]] = []

        for item in path:
            item_name: str = item.identity.matched_rule.name if item.identified else "X"
            item_tags: set[int] = get_tags_mol(item.mol)
            primary_sequence.append((item_name, item_tags))
            item_smiles = mol_to_smiles(item.mol, include_tags=True)
            building_blocks.append(item_smiles)

        product = smiles_to_mol(".".join([x for x in building_blocks]))
        product_smi = mol_to_smiles(product, include_tags=True)
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
