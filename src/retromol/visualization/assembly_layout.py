"""Shared RDKit coordinates for comparing assembly projections of one molecule."""

from dataclasses import dataclass, field
from typing import Iterable

from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.rdchem import Mol


@dataclass(frozen=True)
class MoleculeLayout:
    """A 2D depiction and its original-heavy-atom coordinates, indexed by tag."""

    atom_positions: dict[int, tuple[float, float]]
    bounds: tuple[float, float, float, float]
    molecule: Mol = field(repr=False, compare=False)

    @classmethod
    def from_mol(cls, root_mol: Mol) -> "MoleculeLayout":
        # Restore tag order so serialization/reloading or RDKit atom renumbering
        # cannot rotate/rearrange equivalent projections of the same tagged root.
        order = sorted(range(root_mol.GetNumAtoms()),
                       key=lambda i: (root_mol.GetAtomWithIdx(i).GetIsotope() or float("inf"), i))
        molecule = Chem.RenumberAtoms(root_mol, order) if order else Mol(root_mol)
        tags = [atom.GetIsotope() for atom in molecule.GetAtoms()]
        # Isotopes are RetroMol's atom identifiers, not labels for the depiction.
        for atom in molecule.GetAtoms():
            atom.SetIsotope(0)
            atom.SetAtomMapNum(0)
        positions = {}
        heavy_positions = []
        if molecule.GetNumAtoms():
            rdDepictor.Compute2DCoords(molecule, canonOrient=True, clearConfs=True, forceRDKit=True)
            conformer = molecule.GetConformer()
            for index, atom in enumerate(molecule.GetAtoms()):
                if atom.GetAtomicNum() <= 1:
                    continue
                point = conformer.GetAtomPosition(index)
                xy = (float(point.x), float(point.y))
                heavy_positions.append(xy)
                tag = tags[index]
                if tag:
                    if tag in positions:
                        raise ValueError(f"Root contains duplicate atom tag {tag}")
                    positions[tag] = xy
        if heavy_positions:
            xs, ys = zip(*heavy_positions)
            padding = max(2.0, .12 * max(max(xs) - min(xs), max(ys) - min(ys)))
            bounds = (min(xs) - padding, max(xs) + padding, min(ys) - padding, max(ys) + padding)
        else:
            bounds = (-2.0, 2.0, -2.0, 2.0)
        return cls(positions, bounds, molecule)

    def centroid(self, tags: Iterable[int]) -> tuple[float, float] | None:
        """Centroid of unique original heavy atoms, excluding introduced atoms/H."""
        points = [self.atom_positions[tag] for tag in sorted(set(tags)) if tag in self.atom_positions]
        if not points:
            return None
        return (sum(x for x, _ in points) / len(points), sum(y for _, y in points) / len(points))
