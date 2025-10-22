"""Module for RetroMol results drawing."""

import logging
import os
from copy import deepcopy
from enum import Enum

from rdkit.Chem.Draw.rdMolDraw2D import MolDraw2DSVG, MolDrawOptions

from retromol import chem, config, io, readout


class Palette(Enum):
    Red = (230, 25, 75)
    Blue = (0, 130, 200)
    Green = (60, 180, 75)
    Maroon = (128, 0, 0)
    Brown = (170, 110, 40)
    Olive = (128, 128, 0)
    Teal = (0, 128, 128)
    Navy = (0, 0, 128)
    Orange = (245, 130, 48)
    Yellow = (255, 225, 25)
    Lime = (210, 245, 60)
    Cyan = (70, 240, 240)
    Purple = (145, 30, 180)
    Magenta = (240, 50, 230)
    Pink = (255, 190, 212)
    Apricot = (255, 215, 180)
    Beige = (255, 250, 200)
    Mint = (170, 255, 195)
    Lavender = (220, 190, 255)

    def hex(self, alpha: float) -> str:
        """
        Get hex representation of the color with specified alpha transparency.

        :param alpha: alpha transparency (0.0 to 1.0)
        :return: hex color string with alpha
        """
        return f"#{self.value[0]:02x}{self.value[1]:02x}{self.value[2]:02x}{int(alpha * 255):02x}"

    def normalize(self, min_val: float = 0.0, max_val: float = 255.0) -> tuple[float, float, float]:
        """
        Get normalized RGB tuple of the color.

        :param min_val: minimum value for normalization
        :param max_val: maximum value for normalization
        :return: normalized RGB tuple
        """
        r, g, b = self.value
        return (
            (r - min_val) / (max_val - min_val),
            (g - min_val) / (max_val - min_val),
            (b - min_val) / (max_val - min_val),
        )


def hex_to_rgb_tuple(hex_str: str) -> tuple[float, float, float]:
    """
    Convert hex color string to normalized RGB tuple.

    :param hex_str: hex color string (e.g. "#ff5733" or "#ff5733ff")
    :return: normalized RGB tuple
    """
    hex_str = hex_str.lstrip("#")
    if len(hex_str) == 6:
        r, g, b = int(hex_str[0:2], 16), int(hex_str[2:4], 16), int(hex_str[4:6], 16)
    elif len(hex_str) == 8:
        r, g, b = int(hex_str[0:2], 16), int(hex_str[2:4], 16), int(hex_str[4:6], 16)
        # alpha = int(hex_str[6:8], 16)  # Alpha is ignored in this function
    else:
        raise ValueError(f"Invalid hex color string: {hex_str}")
    return (r / 255.0, g / 255.0, b / 255.0)


def draw_result(
    result: io.Result,
    out_dir: str,
    base_name: str = "optimal_mapping",
    window_size: tuple[int, int] = (800, 800),
    background_color: str | None = None,
) -> None:
    """
    Draw optimal mappings for a RetroMol result.

    For each optimal mapping of identified nodes to the input molecule, generate
    a 2D drawing highlighting the atoms and bonds covered by each node in a
    distinct color.

    :param result: RetroMol Result object containing the input molecule and identified nodes
    :param out_dir: directory to save the output SVG files
    :param base_name: base name for the output files; an index will be appended for each mapping
    :param window_size: size of the drawing window (width, height)
    :param background_color: optional background color for the drawing in
        hex format (e.g. "#ffffff" for white). If None, defaults to transparent
    :return: None
    """
    logger = logging.getLogger(config.LOGGER_NAME)

    optimal_mappings = readout.optimal_mappings_with_timeout(result)
    logger.info(f"{len(optimal_mappings)} optimal mapping(s) found.")

    # Retrieve input SMILES with tags
    input_smi = result.get_input_smiles(remove_tags=False)
    input_mol = chem.smiles_to_mol(input_smi)

    # If no optimal mappings, just draw the input molecule
    if len(optimal_mappings) == 0:
        optimal_mappings.append({})

    for map_idx, o_m in enumerate(optimal_mappings):
        out_path = os.path.join(out_dir, f"{base_name}_{map_idx + 1}.svg")

        drawing: MolDraw2DSVG = MolDraw2DSVG(*window_size)
        palette = [c.normalize() for c in Palette]

        atoms_to_highlight: list[int] = []
        bonds_to_highlight: list[int] = []
        atom_highlight_colors: dict[int, tuple[float, float, float]] = {}
        bond_highlight_colors: dict[int, tuple[float, float, float]] = {}

        for n_idx, node in enumerate(o_m.get("nodes", [])):
            color = palette[n_idx % len(palette)]
            n_tags = node["tags"]

            logger.info(f"Mapping {map_idx + 1} - node {n_idx + 1}: {node['identity']} {n_tags}")

            for atom in input_mol.GetAtoms():
                a_tag = atom.GetIsotope()
                if a_tag in n_tags:
                    a_idx = atom.GetIdx()
                    atoms_to_highlight.append(a_idx)
                    atom_highlight_colors[a_idx] = color

            for bond in input_mol.GetBonds():
                b_begin_idx = bond.GetBeginAtom()
                b_end_idx = bond.GetEndAtom()
                b_begin_tag = b_begin_idx.GetIsotope()
                b_end_tag = b_end_idx.GetIsotope()
                if b_begin_tag in n_tags and b_end_tag in n_tags:
                    b_idx = bond.GetIdx()
                    bonds_to_highlight.append(b_idx)
                    bond_highlight_colors[b_idx] = color

        options: MolDrawOptions = drawing.drawOptions()
        if background_color is not None:
            options.setBackgroundColour(hex_to_rgb_tuple(background_color))
        options.useBWAtomPalette()

        # Remove isotopic labels for drawing
        cp_input_mol = deepcopy(input_mol)
        for atom in cp_input_mol.GetAtoms():
            atom.SetIsotope(0)

        drawing.DrawMolecule(
            cp_input_mol,
            highlightAtoms=atoms_to_highlight,
            highlightBonds=bonds_to_highlight,
            highlightAtomColors=atom_highlight_colors,
            highlightBondColors=bond_highlight_colors,
        )

        drawing.FinishDrawing()
        svg_str = drawing.GetDrawingText().replace("svg:", "")

        with open(out_path, "w") as f:
            f.write(svg_str)
        logger.info(f"Wrote drawing to {out_path}")

        # Close the drawing to free up memory
        drawing = None
