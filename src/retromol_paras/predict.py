"""Runs the trained PARAS model against new adenylation-domain sequences."""

from dataclasses import dataclass
from pathlib import Path
from tempfile import TemporaryDirectory

from sklearn.ensemble import RandomForestClassifier

from retromol_paras.constants import SMILES_FILE
from retromol_paras.domain import AdenylationDomain
from retromol_paras.featurisation import extract_domains, get_domain_features


@dataclass(frozen=True)
class SubstratePrediction:
    """ One (label, probability) call for a single adenylation domain. """

    label: str
    score: float
    smiles: str | None


@dataclass(frozen=True)
class DomainPrediction:
    """ All substrate predictions for one detected adenylation domain, best first. """

    protein_name: str
    domain_nr: int
    start: int
    end: int
    extended_signature: str
    predictions: list[SubstratePrediction]


def _load_label_to_smiles() -> dict[str, str]:
    label_to_smiles: dict[str, str] = {}
    with open(SMILES_FILE, "r") as fo:
        next(fo)  # header
        for line in fo:
            name, smiles = line.rstrip("\n").split("\t")
            label_to_smiles[name] = smiles
    return label_to_smiles


_LABEL_TO_SMILES: dict[str, str] | None = None


def label_to_smiles(label: str) -> str | None:
    """ Look up a predicted substrate label's SMILES from PARAS' packaged substrate table. """
    global _LABEL_TO_SMILES
    if _LABEL_TO_SMILES is None:
        _LABEL_TO_SMILES = _load_label_to_smiles()
    return _LABEL_TO_SMILES.get(label)


def predict_domains(
    model: RandomForestClassifier,
    protein_sequences: dict[str, str],
    keep_top: int = 3,
    path_temp_dir: str | Path | None = None,
) -> list[DomainPrediction]:
    """
    Detect and predict the substrate of every adenylation domain in a set of protein sequences.

    :param model: a fitted PARAS RandomForestClassifier (see retromol_paras.train.train_model).
    :param protein_sequences: {sequence_id: amino_acid_sequence}.
    :param keep_top: number of top-scoring substrate predictions to keep per domain.
    :param path_temp_dir: scratch directory for HMMER/MUSCLE intermediates; a temp dir is used if not given.
    :return: one DomainPrediction per detected adenylation domain.
    """
    def _run(scratch: Path) -> list[DomainPrediction]:
        fasta_file = scratch / "query.fasta"
        with open(fasta_file, "w") as fo:
            for seq_id, seq in protein_sequences.items():
                fo.write(f">{seq_id}\n{seq}\n")

        domains: list[AdenylationDomain] = extract_domains(fasta_file, scratch)
        if not domains:
            return []

        features = [get_domain_features(d.extended_signature) for d in domains]
        probabilities = model.predict_proba(features)
        class_labels = model.classes_.tolist()

        results = []
        for domain, probs in zip(domains, probabilities):
            ranked = sorted(zip(class_labels, probs), key=lambda x: x[1], reverse=True)[:keep_top]
            results.append(DomainPrediction(
                protein_name=domain.protein_name,
                domain_nr=domain.domain_nr,
                start=domain.start,
                end=domain.end,
                extended_signature=domain.extended_signature,
                predictions=[
                    SubstratePrediction(label=label, score=float(prob), smiles=label_to_smiles(label))
                    for label, prob in ranked
                ],
            ))
        return results

    if path_temp_dir is not None:
        return _run(Path(path_temp_dir))

    with TemporaryDirectory(prefix="retromol_paras_") as tmp:
        return _run(Path(tmp))
