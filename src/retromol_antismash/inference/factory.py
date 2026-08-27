"""
Single place that turns pmp.yml's predictor selection into an actual, ready-to-register
DomainInferenceModel -- so the database pipeline (database/scripts/parse_gbks.py), the
GUI backend (gui/src/server/routes/jobs.py), and the CLI test script
(scripts/predict_bgc.py) all pick their NRPS A-domain substrate predictor the same way,
from the same config, instead of each hardcoding its own model choice. Switching
predictors (or adding a new one) is a pmp.yml edit plus one branch here, not a
multi-file change.

There is deliberately no pyhmmer-based model here: "paras_cli"
(retromol_antismash.inference.model_paras_cli.ParasCliModel, backed by
src/retromol_paras/) is the only NRPS substrate-prediction model, running domain
extraction via command-line HMMER2/HMMER3/MUSCLE3 and training its own
RandomForestClassifier against this repo's pinned scikit-learn version, rather than
unpickling a pretrained model file built under someone else's environment.
"""

from pathlib import Path

from retromol_antismash.inference.base import DomainInferenceModel
from retromol_antismash.predictions import PredictionConfig


def build_nrps_a_domain_model(
    config: PredictionConfig | None = None,
    *,
    threshold: float = 0.1,
    keep_top: int = 3,
    cache_dir: str | Path | None = None,
    force_retrain: bool = False,
) -> DomainInferenceModel | None:
    """
    Build (but don't register -- see
    retromol_antismash.inference.registry.register_domain_model) whichever model
    pmp.yml's "predictors.nrps.a_domain" currently selects.

    :param config: prediction configuration to read the selection from, or None to
        load the packaged default pmp.yml.
    :param threshold: minimum predicted probability for a substrate call to be kept
        (passed through to whichever model is built).
    :param keep_top: number of top-scoring substrate predictions to keep per domain.
    :param cache_dir: training-signature + fitted-model cache directory (see
        retromol_paras.train.train_model).
    :param force_retrain: retrain from scratch even if a cached model is present.
    :return: a DomainInferenceModel, or None if the selected method reads
        antiSMASH's own GenBank qualifier instead of running a model
        (source: qualifier) -- there's nothing to build or register in that case;
        retromol_antismash.modules.collect_nrps_modules reads the qualifier
        directly off the domain when it resolves the module's substrate.
    :raises ValueError: if the selected method is source: model but names a model
        this factory doesn't know how to build.
    """
    config = config or PredictionConfig.load_default()
    method = config.get_method("nrps.a_domain")

    if method.source != "model":
        return None

    if method.model == "paras_cli":
        from retromol_antismash.inference.model_paras_cli import ParasCliModel
        return ParasCliModel(threshold=threshold, keep_top=keep_top, cache_dir=cache_dir, force_retrain=force_retrain)

    raise ValueError(
        f"pmp.yml's predictors.nrps.a_domain selects model {method.model!r}, but no "
        "factory is registered for it (see retromol_antismash.inference.factory)"
    )
