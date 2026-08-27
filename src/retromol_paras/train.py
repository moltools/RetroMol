"""
Trains (and caches) the PARAS random-forest classifier from scratch, against
`retromol_paras`'s own scikit-learn version rather than trusting a pickled
model file trained under someone else's environment.

Reproduces `parasect.model_training.rf.train_rf.train_paras_signatures` /
`train_random_forest`, but reading training examples from the flat
`parasect_dataset.txt` file (see constants.PARASECT_DATASET_FILE) instead of
PARASECT's SQLAlchemy database -- no database needed.
"""

import logging
from pathlib import Path

import joblib
import numpy as np
from sklearn.ensemble import RandomForestClassifier

from retromol_paras.constants import PARASECT_DATASET_FILE
from retromol_paras.featurisation import extract_domains, get_domain_features
from retromol_paras.tabular import Tabular

log = logging.getLogger(__name__)

RANDOM_STATE = 25051989  # matches parasect.model_training.rf.train_rf.train_random_forest
N_ESTIMATORS = 1000

DEFAULT_CACHE_DIR = Path.home() / ".retromol_paras"
SIGNATURE_CACHE_FILE = "extended_signatures.cache.tsv"
MODEL_FILE = "model.paras.joblib"


def _parse_dataset(path_in: Path) -> tuple[dict[str, str], dict[str, str]]:
    """ Read parasect_dataset.txt into {domain_id: sequence} and {domain_id: specificity}. """
    data = Tabular(path_in, separator="\t")
    id_to_seq = {row_id: data.get_row_value(row_id, "sequence") for row_id in data.rows}
    # a domain_id can carry several pipe-separated synonyms mapping to the same specificity;
    # parasect just uses the whole id string as the training label's key, so we do too.
    id_to_spec = {row_id: data.get_row_value(row_id, "specificity") for row_id in data.rows}
    return id_to_seq, id_to_spec


def _extract_training_signatures(cache_dir: Path, force: bool = False) -> dict[str, str]:
    """
    Run every parasect_dataset.txt domain through the HMMER2/HMMER3/MUSCLE3 extraction
    pipeline once, caching {domain_id: extended_signature} to `cache_dir` since this is
    the slow step (thousands of subprocess calls to legacy tools) and doesn't change
    unless the packaged dataset does.
    """
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_file = cache_dir / SIGNATURE_CACHE_FILE

    if cache_file.exists() and not force:
        log.info(f"reusing cached extended signatures: {cache_file}")
        id_to_signature: dict[str, str] = {}
        with open(cache_file, "r") as fo:
            next(fo)  # header
            for line in fo:
                domain_id, signature = line.rstrip("\n").split("\t")
                id_to_signature[domain_id] = signature
        return id_to_signature

    id_to_seq, _ = _parse_dataset(PARASECT_DATASET_FILE)

    scratch = cache_dir / "extraction_scratch"
    scratch.mkdir(parents=True, exist_ok=True)
    dataset_fasta = scratch / "dataset.fasta"
    with open(dataset_fasta, "w") as fo:
        for domain_id, seq in id_to_seq.items():
            fo.write(f">{domain_id}\n{seq}\n")

    log.info(f"extracting extended signatures for {len(id_to_seq)} training domains (this runs HMMER2/HMMER3/MUSCLE3 once, then caches) ...")
    domains = extract_domains(dataset_fasta, scratch)

    # one label per training row (domain_id); keep the first detected A-domain per header.
    id_to_signature = {}
    for domain in domains:
        if domain.protein_name not in id_to_signature and domain.extended_signature:
            id_to_signature[domain.protein_name] = domain.extended_signature

    with open(cache_file, "w") as fo:
        fo.write("domain_id\textended_signature\n")
        for domain_id, signature in sorted(id_to_signature.items()):
            fo.write(f"{domain_id}\t{signature}\n")

    return id_to_signature


def train_model(cache_dir: str | Path | None = None, force: bool = False) -> RandomForestClassifier:
    """
    Train (or load a cached) PARAS random-forest classifier.

    :param cache_dir: directory to cache the extracted training signatures and the fitted
        model in. Defaults to ~/.retromol_paras.
    :param force: retrain from scratch even if a cached model is present.
    :return: a fitted RandomForestClassifier, predicting substrate name from a 34-residue
        extended signature's physicochemical feature vector (see featurisation.get_domain_features).
    """
    cache_dir = Path(cache_dir) if cache_dir else DEFAULT_CACHE_DIR
    cache_dir.mkdir(parents=True, exist_ok=True)
    model_file = cache_dir / MODEL_FILE

    if model_file.exists() and not force:
        log.info(f"loading cached model: {model_file}")
        return joblib.load(model_file)

    id_to_seq, id_to_spec = _parse_dataset(PARASECT_DATASET_FILE)
    id_to_signature = _extract_training_signatures(cache_dir, force=force)

    features: list[list[float]] = []
    labels: list[str] = []
    for domain_id in id_to_seq:
        signature = id_to_signature.get(domain_id)
        if not signature:
            continue
        features.append(get_domain_features(signature))
        # parasect_dataset.txt can list multiple pipe-separated specificities for a promiscuous
        # domain; PARAS' default (FIRST_ONLY) selection mode trains on just the first.
        labels.append(id_to_spec[domain_id].split("|")[0])

    log.info(f"training PARAS classifier on {len(labels)} domains ...")
    model = RandomForestClassifier(
        n_estimators=N_ESTIMATORS,
        n_jobs=1,
        oob_score=True,
        random_state=RANDOM_STATE,
        class_weight="balanced",
    )
    model.fit(np.array(features), np.array(labels))

    joblib.dump(model, model_file)
    log.info(f"cached trained model: {model_file} (oob_score={model.oob_score_:.4f})")

    return model
