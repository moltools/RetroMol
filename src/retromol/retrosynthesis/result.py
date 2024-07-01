# -*- coding: utf-8 -*-

"""This module contains the dataclass for storing the results of a RetroMol run."""

import json
import typing as ty
from dataclasses import dataclass

from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Geometry import Point2D

from retromol.retrosynthesis.drawing import Palette, get_2d_coordinatates


@dataclass
class Result:
    """Dataclass for storing the results of a RetroMol run.

    :param identifier: The identifier of the molecule.
    :type identifier: str
    :param mol: The molecule.
    :type mol: Chem.Mol
    :param success: Whether the RetroMol run was successful.
    :type success: bool
    :param reaction_tree: The reaction tree.
    :type reaction_tree: ty.Dict[int, ty.List[int]]
    :param applied_reactions: The applied reactions.
    :type applied_reactions: ty.List[str]
    :param monomer_graph: The monomer graph.
    :type monomer_graph: ty.Dict[str, ty.Any]
    :param sequences: The modular natural product sequences.
    :type sequences: ty.List[ty.List[str]]
    """

    identifier: str
    mol: Chem.Mol
    success: bool
    reaction_tree: ty.Dict[int, ty.List[int]] = None
    applied_reactions: ty.List[str] = None
    monomer_graph: ty.Dict[str, ty.Any] = None
    sequences: ty.List[ty.List[str]] = None

    def to_json(self, indent: int = 4) -> str:
        """Convert the result to JSON.

        :param indent: The number of spaces to indent the JSON.
        :type indent: int
        :return: The result as JSON.
        :rtype: str
        """
        return json.dumps(
            {
                "identifier": self.identifier,
                "smiles": Chem.MolToSmiles(self.mol),
                "success": "true" if self.success else "false",
                "reaction_tree": self.reaction_tree,
                "applied_reactions": self.applied_reactions,
                "monomer_graph": self.monomer_graph,
                "sequences": self.sequences,
            },
            indent=indent,
        )

    @classmethod
    def from_json(cls, path: str) -> "Result":
        """Create a Result object from a JSON file.

        :param path: The path to the JSON file.
        :type path: str
        :return: The Result object.
        :rtype: Result
        """
        with open(path, "r", encoding="utf-8") as fo:
            data = json.load(fo)

        reaction_tree = data["reaction_tree"]
        if reaction_tree is not None:
            reaction_tree = {int(k): v for k, v in reaction_tree.items()}

        monomer_graph = data["monomer_graph"]
        if monomer_graph is not None:
            monomer_graph = {int(k): v for k, v in monomer_graph.items()}

        return Result(
            identifier=data["identifier"],
            mol=Chem.MolFromSmiles(data["smiles"]),
            success=True if data["success"] == "true" else False,
            reaction_tree=reaction_tree,
            applied_reactions=data["applied_reactions"],
            monomer_graph=monomer_graph,
            sequences=data["sequences"],
        )

    def has_identified_monomers(self) -> bool:
        """Check if any monomers have been identified.

        :return: Whether any monomers have been identified.
        :rtype: bool
        """
        return any([x["identity"] is not None for _, x in self.monomer_graph.items()])
    
    def draw_molecule(
        self, 
        color_monomers: bool = False,
        window_size: ty.Tuple[int, int] = (800, 800),
        background_color: ty.Optional[str] = None,
    ) -> str:
        """Draw the monomer graph.

        :param color_monomers: Whether to color the monomers in the molecule.
        :type color_monomers: bool
        :param window_size: The window size.
        :type window_size: ty.Tuple[int, int]
        :param background_color: The background color.
        :type background_color: ty.Optional[str]
        :return: The SVG string.
        :rtype: str
        """
        # Reconstruct mol.
        mol = self.mol

        # Create atom map number to atom index mapping.
        mapping = {}
        for atom in mol.GetAtoms():
            atom_map_num = atom.GetAtomMapNum()
            if atom_map_num != 0:
                mapping[atom_map_num] = atom.GetIdx()

        # Color identifier monomers in the molecule.
        drawing = rdMolDraw2D.MolDraw2DSVG(*window_size)
        palette = [c.normalize() for c in Palette]

        atoms_to_highlight = []
        bonds_to_highlight = []
        atom_highlight_colors = {}
        bond_highlight_colors = {}
        labels = []

        if color_monomers:
            coords = get_2d_coordinatates(mol)
            monomer_count = 0

            # Get all monomer atom mappings and their identity.
            for _, monomer in self.monomer_graph.items():

                if monomer_identity := monomer["identity"]:
                    atom_indices = []
                    color = palette[monomer_count % len(palette)]

                    # Get atom indices.
                    reaction_tree_id = monomer["reaction_tree_id"]
                    monomer_smiles = self.reaction_tree[reaction_tree_id]["smiles"]
                    monomer_mol = Chem.MolFromSmiles(monomer_smiles)
                    monomer_inds = [
                        mapping[x.GetAtomMapNum()]
                        for x in monomer_mol.GetAtoms() 
                        if x.GetAtomMapNum() != 0
                    ]

                    # Get centroid for monomer.
                    x_coords = [coords[x][0] for x in monomer_inds]
                    y_coords = [coords[x][1] for x in monomer_inds]
                    centroid = (sum(x_coords) / len(x_coords), sum(y_coords) / len(y_coords))
                    centroid = Point2D(*centroid)
                    labels.append((monomer_identity, centroid))
                    
                    # Highlight atoms.
                    for atom_idx in monomer_inds:
                        atom_indices.append(atom_idx)
                        atoms_to_highlight.append(atom_idx)
                        atom_highlight_colors[atom_idx] = color

                    # Highlight bonds.
                    for bond in mol.GetBonds():
                        bond_index = bond.GetIdx()
                        start_atom = bond.GetBeginAtom().GetIdx()
                        end_atom = bond.GetEndAtom().GetIdx()

                        if start_atom in atom_indices and end_atom in atom_indices:
                            bonds_to_highlight.append(bond_index)
                            bond_highlight_colors[bond_index] = color

                    monomer_count += 1

        options = drawing.drawOptions()
        if background_color is not None:
            options.setBackgroundColour(background_color)
        options.useBWAtomPalette()

        # Make sure atom map numbers are not drawn.
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(0)

        drawing.DrawMolecule(
            mol,
            highlightAtoms=atoms_to_highlight,
            highlightBonds=bonds_to_highlight,
            highlightAtomColors=atom_highlight_colors,
            highlightBondColors=bond_highlight_colors
        )

        # Draw labels.
        for label, position in labels:
            drawing.DrawString(label, position)

        drawing.FinishDrawing()
        svg_str = drawing.GetDrawingText().replace("svg:", "")

        return svg_str
