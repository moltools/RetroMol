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


rxn_start_capping = smarts_to_reaction("[SnH:1][C:2]~[C:3][PbH:4]>>[OH][C](=O)[C:2]~[C:3][PbH:4].[SnH2:1]")
rxn_extension_1 = smarts_to_reaction("[C:1][C:2]~[C:3][PbH:4].[SnH:5][C:6]-[C:7][PbH:8]>>[C:1][C:2]~[C:3][C:6]-[C:7][PbH:8].[PbH2:4].[SnH2:5]")
rxn_extension_2 = smarts_to_reaction("[C:1][C:2]~[C:3][PbH:4].[SnH:5][C:6]=[C:7][PbH:8]>>[C:1][C:2]~[C:3][C:6]=[C:7][PbH:8].[PbH2:4].[SnH2:5]")
# only works for existing polyketide chain, extender amino acid
# also need:
# - existing amino acid chain, extender amino acid
# - existing amino acid chain, extender polyketide unit
rxn_extension_3 = smarts_to_reaction("[C:1][C:2]~[C:3][PbH:4].[OH:5][C:6](=[O:7])[C:8][NH2:9]>>[C:1][C:2]~[C:3][C:8][NH2:9].[PbH2:4].[O:5]=[C:6](=[O:7])")
extension_1_pattern = smarts_to_mol("[SnH][C]-[C][PbH]")
extension_2_pattern = smarts_to_mol("[SnH][C]=[C][PbH]")
extension_3_pattern = smarts_to_mol("[OH][C](=O)C[NH2]")
rxn_end_capping = smarts_to_reaction("[PbH:1][C:2]>>[*][C:2].[PbH2:1]")


indicator_smarts = ["[PbH2]", "[SnH2]", "[O]=[C]=[O]"]
indicator_mols = [smarts_to_mol(s) for s in indicator_smarts]
def filter_indicators(products):
    results = []
    for product in products:
        if not any(product.HasSubstructMatch(ind) for ind in indicator_mols):
            results.append(product)
    return results


def reconstruct_linear_readout(result: Result) -> Reconstruction:
    """
    Reconstruct a linear backbone from RetroMol's linear readout.

    :param result: The result object returned by RetroMol, containing the reaction graph and the original molecule.
    :return: The reconstructed backbone.
    """
    root_enc = encode_mol(result.submission.mol)
    readout = LinearReadout.from_reaction_graph(root_enc, reaction_graph=result.reaction_graph, identified_only=True)

    # Get building blocks
    paths = readout.paths
    paths.sort(key=lambda x: len(x), reverse=True)
    building_blocks = []
    primary_sequence: list[tuple[str, set[int]]] = []
    for path in readout.paths:
        for item in path:
            item_name: str = item.identity.matched_rule.name if item.identified else "X"
            item_tags: set[int] = get_tags_mol(item.mol)
            primary_sequence.append((item_name, item_tags))
            item_smiles = mol_to_smiles(item.mol, include_tags=True)
            building_block = smiles_to_mol(item_smiles)
            building_blocks.append(building_block)
        break

    # Start reconstruction
    # TODO: capping start could also be an amino acid!
    products = rxn_start_capping.RunReactants((building_blocks[0],))
    assert len(products) == 1, "Expected exactly one product!"
    products = list(products[0])
    products = filter_indicators(products)
    assert len(products) == 1, "Expected exactly one product!"
    product = products[0]

    # Loop through remaining items and attach to first building block
    for building_block in building_blocks[1:]:
        print(get_tags_mol(building_block))
        if building_block.HasSubstructMatch(extension_1_pattern):
            products = rxn_extension_1.RunReactants((product, building_block,))
        elif building_block.HasSubstructMatch(extension_2_pattern):
            products = rxn_extension_2.RunReactants((product, building_block,))
        elif building_block.HasSubstructMatch(extension_3_pattern):
            products = rxn_extension_3.RunReactants((product, building_block,))
        else:
            raise ValueError(f"Unknown building block: {building_block}")

        products = products[0]
        products = filter_indicators(products)
        assert len(products) == 1, "Expected exactly one product!"
        product = products[0]

    # End reconstruction; cap
    try:
        products = rxn_end_capping.RunReactants((product,))
        assert len(products) == 1, "Expected exactly one product!"
        products = filter_indicators(products[0])
        product = products[0]
    except:
        pass

    product_smi = mol_to_smiles(product, include_tags=True)
    print(product_smi)

    input_smi = mol_to_smiles(result.submission.mol, include_tags=True)
    print(input_smi)

    print(primary_sequence)

    return Reconstruction(
        tagged_input_smiles=input_smi,
        tagged_backbone_smiles=product_smi,
        primary_sequence=primary_sequence,
    )
