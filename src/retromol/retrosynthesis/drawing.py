# -*- coding: utf-8 -*-

"""This module contains functions for drawing."""

from enum import Enum
import typing as ty

from rdkit import Chem
from rdkit.Chem import AllChem


class Palette(Enum):
    """Palette of colors for drawing molecules as RGB.
    """
    Red         = (230,  25,  75); Blue        = (  0, 130, 200); Green       = ( 60, 180,  75); Maroon      = (128,   0,   0)
    Brown       = (170, 110,  40); Olive       = (128, 128,   0); Teal        = (  0, 128, 128); Navy        = (  0,   0, 128)
    Orange      = (245, 130,  48); Yellow      = (255, 225,  25); Lime        = (210, 245,  60); Cyan        = ( 70, 240, 240)
    Purple      = (145,  30, 180); Magenta     = (240,  50, 230); Pink        = (255, 190, 212); Apricot     = (255, 215, 180)
    Beige       = (255, 250, 200); Mint        = (170, 255, 195); Lavender    = (220, 190, 255)

    def as_hex(self, alpha: ty.Optional[float] = None) -> str:
        """Return the hex code of the color with the given alpha value.

        :param alpha: The alpha value of the color.
        :type alpha: ty.Optional[float]
        :return: The hex code of the color.
        :rtype: str
        """
        r, g, b = self.value
        hex_base = "#{:02x}{:02x}{:02x}".format(r, g, b)
        if alpha is not None:
            if alpha < 0: alpha = 0.0
            elif alpha > 1: alpha = 1.0
            return "{}{:02x}".format(hex_base, int(alpha * 255))
        else:
            return hex_base

    def normalize(self, minv: float = 0.0, maxv: float  = 255.0) -> ty.Tuple[float, float, float]:
        """Normalize the color values to the range [0, 1].

        :param minv: The minimum value.
        :type minv: float
        :param maxv: The maximum value.
        :type maxv: float
        :return: The normalized color values.
        :rtype: ty.Tuple[float, float, float]
        """
        r, g, b = self.value
        return (
            round((r-minv)/(maxv-minv), 3), 
            round((g-minv)/(maxv-minv), 3), 
            round((b-minv)/(maxv-minv), 3)
        )


def get_2d_coordinatates(mol: Chem.Mol) -> ty.Dict[int, ty.Tuple[float, float]]:
    """Get 2D coordinates of atoms in molecule.

    :param mol: The molecule to get the 2D coordinates of.
    :type mol: Chem.Mol
    :return: A dictionary of atom tracking numbers and their 2D coordinates.
    :rtype: ty.Dict[int, ty.Tuple[float, float]]
    :note: The atom tracking numbers are stored as atom isotope numbers.
    """
    AllChem.Compute2DCoords(mol)

    coordinates = {}
    for atom in mol.GetAtoms():
        atom_idx = atom.GetIdx()
        position = mol.GetConformer(0).GetAtomPosition(atom_idx)
        coordinates[atom_idx] = (position.x, position.y)

    return coordinates
