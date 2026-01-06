"""Module defining the Submission data class."""

from dataclasses import dataclass, field
from typing import Any

from rdkit.Chem.rdchem import Mol

from retromol.chem.mol import standardize_from_smiles, mol_to_inchikey
from retromol.chem.tagging import tag_mol


@dataclass(frozen=True)
class Submission:
    """
    Represents a submission of a molecule for retrosynthetic analysis.

    :var smiles: str: SMILES representation of the submitted molecule
    :var name: str | None: optional name of the submitted molecule
    :var props: dict[str, Any] | None: optional additional properties associated with the submission
    :var keep_stereo: bool: whether to keep stereochemistry during standardization
    :var neutralize: bool: whether to neutralize the molecule during standardization
    :var canonicalize_tautomer: bool: whether to canonicalize the tautomer during
    """

    smiles: str
    name: str | None = None
    props: dict[str, Any] | None = None

    keep_stereo: bool = True
    neutralize: bool = True
    canonicalize_tautomer: bool = False

    mol: Mol = field(init=False, repr=False)
    inchikey: str = field(init=False, repr=False)

    def __post_init__(self) -> None:
        """
        Post-initialization processing to generate standardized molecule and InChIKey.
        """
        # Sanitize SMILES
        smiles = self.smiles.replace("[N]", "N")  # avoid parsing issues with RDKit

        # Generate standardized molecule
        mol = standardize_from_smiles(
            smiles,
            keep_stereo=self.keep_stereo,
            neutralize=self.neutralize,
            tautomer_canon=self.canonicalize_tautomer,
        )

        # Generate InChIKey
        inchikey = mol_to_inchikey(mol)

        # Tag molecule
        tag_mol(mol)

        object.__setattr__(self, "smiles", smiles)
        object.__setattr__(self, "mol", mol)
        object.__setattr__(self, "inchikey", inchikey)

    def __str__(self) -> str:
        """
        String representation of the Submission.

        :return: string representation of the Submission
        """
        return f"Submission(name={self.name})"
    
    def to_dict(self) -> dict[str, Any]:
        """
        Serialize the Submission to a dictionary.

        :return: dictionary representation of the Submission
        """
        return {
            "smiles": self.smiles,
            "name": self.name,
            "props": self.props,
        }
    
    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "Submission":
        """
        Deserialize a Submission from a dictionary.

        :param data: dictionary representation of the Submission
        :return: Submission object
        """
        return cls(
            smiles=data["smiles"],
            name=data.get("name"),
            props=data.get("props"),
        )
