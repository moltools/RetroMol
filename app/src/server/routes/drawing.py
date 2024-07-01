"""This module contains the API endpoints for CineMol."""
from flask import Blueprint, Response, request
from enum import Enum
import typing as ty

from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

from .common import ResponseStatus, ResponseData


class Palette(Enum):
    Red         = (230,  25,  75)
    Blue        = (  0, 130, 200)
    Green       = ( 60, 180,  75)
    Maroon      = (128,   0,   0)
    Brown       = (170, 110,  40)
    Olive       = (128, 128,   0)
    Teal        = (  0, 128, 128)
    Navy        = (  0,   0, 128)
    Orange      = (245, 130,  48)
    Yellow      = (255, 225,  25)
    Lime        = (210, 245,  60)
    Cyan        = ( 70, 240, 240)
    Purple      = (145,  30, 180)
    Magenta     = (240,  50, 230)
    Pink        = (255, 190, 212)
    Apricot     = (255, 215, 180)
    Beige       = (255, 250, 200)
    Mint        = (170, 255, 195)
    Lavender    = (220, 190, 255)

    def hex(self, alpha: float) -> str:
        return "#{:02x}{:02x}{:02x}{:02x}".format(
            self.value[0], 
            self.value[1], 
            self.value[2],
            int(alpha * 255)
        )

    def normalize(
        self,
        min_val: ty.Optional[float] = 0.0,
        max_val: ty.Optional[float] = 255.0
    ) -> ty.Tuple[float, float, float]:
        r, g, b = self.value
        return (
            (r - min_val) / (max_val - min_val), 
            (g - min_val) / (max_val - min_val), 
            (b - min_val) / (max_val - min_val)
        )


def draw(
    mol: Chem.Mol,
    subs: ty.List[ty.List[int]] = [],
    window_size: ty.Tuple[int, int] = (500, 500),
    background_color: ty.Optional[str] = None
) -> str:
    """
    Draw a molecule with its substructures highlighted.

    Args:
        mol (Chem.Mol): Molecule.
        subs (ty.List[Chem.Mol], optional): List of substructure indices to highlight. Defaults to [].
        window_size (ty.Tuple[int, int], optional): Window size. Defaults to (800, 800).
        background_color (ty.Optional[str], optional): Background color. Defaults to None.

    Returns:
        str: SVG string.
    """
    drawing = rdMolDraw2D.MolDraw2DSVG(*window_size)
    palette = [c.normalize() for c in Palette]

    amn_to_idx = {}
    for atom in mol.GetAtoms():
        amn = atom.GetAtomMapNum()
        idx = atom.GetIdx()
        if amn > 0:
            amn_to_idx[amn] = idx

    atoms_to_highlight = []
    bonds_to_highlight = []
    atom_highlight_colors = {}
    bond_highlight_colors = {}

    for i, sub in enumerate(subs):
        atom_indices = []
        color = palette[i % len(palette)]

        for amn in sub:
            atom_index = amn_to_idx[amn]
            atom_indices.append(atom_index)
            atoms_to_highlight.append(atom_index)
            atom_highlight_colors[atom_index] = color

        for bond in mol.GetBonds():
            bond_index = bond.GetIdx()
            start_atom = bond.GetBeginAtom().GetIdx()
            end_atom = bond.GetEndAtom().GetIdx()

            if start_atom in atom_indices and end_atom in atom_indices:
                bonds_to_highlight.append(bond_index)
                bond_highlight_colors[bond_index] = color

    options = drawing.drawOptions()
    if background_color is not None:
        options.setBackgroundColour(background_color)
    options.useBWAtomPalette()

    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)

    drawing.DrawMolecule(
        mol,
        highlightAtoms=atoms_to_highlight,
        highlightBonds=bonds_to_highlight,
        highlightAtomColors=atom_highlight_colors,
        highlightBondColors=bond_highlight_colors
    )

    drawing.FinishDrawing()
    svg_str = drawing.GetDrawingText().replace("svg:", "")

    return svg_str


blueprint_draw_smiles_with_highlights = Blueprint("draw_smiles_with_highlights", __name__)
@blueprint_draw_smiles_with_highlights.route("/api/draw_smiles_with_highlights", methods=["POST"])
def draw_smiles_with_highlights() -> Response:
    """API endpoint for drawing SMILES with highlights."""
    # Parse request data.
    data = request.get_json()

    try:
        smiles = data["smiles"]
        highlights = [x[1] for x in data["highlights"]]
        dimensions = data["dimensions"]
        window_size = (dimensions["width"] - 25, dimensions["height"] - 25)

        svg_string = draw(Chem.MolFromSmiles(smiles), highlights, window_size)

        message = "Successfully created SVG with highlights of SMILES!"
        payload =  {"svgString": svg_string}
        return ResponseData(ResponseStatus.Success, payload, message).jsonify()

    except Exception as e:
        msg = f"Failed to draw SVG: {e}"
        return ResponseData(ResponseStatus.Failure, message=msg).jsonify()
    

blueprint_draw_smiles = Blueprint("draw_smiles", __name__)
@blueprint_draw_smiles.route("/api/draw_smiles", methods=["POST"])
def draw_smiles() -> Response:
    """API endpoint for drawing SMILES."""
    data = request.get_json()

    try:
        smiles = data["smiles"]
        dimensions = data["dimensions"]
        window_size = (dimensions["width"], dimensions["height"])
        svg_string = draw(Chem.MolFromSmiles(smiles), [], window_size)
        message = "Successfully created SVG of SMILES!"
        payload = {"svgString": svg_string}
        return ResponseData(ResponseStatus.Success, payload, message).jsonify()
    
    except Exception as e:
        msg = f"Failed to draw SVG: {e}"
        return ResponseData(ResponseStatus.Failure, message=msg).jsonify()
