"""
Predicts NRPS A-domain substrate specificity using retromol_paras, a self-contained
reimplementation of PARAS that runs domain extraction via command-line HMMER2/HMMER3/
MUSCLE3 (no pyhmmer) and trains its own RandomForestClassifier against this repo's
pinned scikit-learn version (no unpickling a pretrained model file built under
someone else's environment).

This is the only NRPS substrate-prediction model in this repo -- registered as
DomainInferenceModel.name "paras_cli", selected by pmp.yml's
"predictors.nrps.a_domain" (see retromol_paras/__init__.py for the required system
binaries, and retromol_antismash.inference.factory for how pmp.yml resolves to
this class).
"""

import logging
from pathlib import Path

from retromol_antismash.inference.base import DomainInferenceModel
from retromol_antismash.model import Domain, InferenceResult
from retromol_paras.predict import predict_domains
from retromol_paras.train import DEFAULT_CACHE_DIR, train_model

log = logging.getLogger(__name__)


class ParasCliModel(DomainInferenceModel):
    """ CLI-tool-based, self-trained reimplementation of the PARAS A-domain substrate predictor. """

    name: str = "paras_cli"

    def __init__(
        self,
        threshold: float = 0.1,
        keep_top: int = 3,
        cache_dir: Path | str | None = None,
        force_retrain: bool = False,
        use_muscle_fallback: bool = True,
    ) -> None:
        """
        :param threshold: minimum predicted probability for a substrate call to be kept.
        :param keep_top: number of top-scoring substrate predictions to keep per domain.
        :param cache_dir: directory to cache the trained model and extracted training
            signatures in (see retromol_paras.train.train_model). Defaults to ~/.retromol_paras.
        :param force_retrain: retrain from scratch even if a cached model is present.
        :param use_muscle_fallback: see retromol_paras.featurisation.get_domains -- False
            trades some coverage (domains HMMER2 misses go unpredicted instead of being
            recovered via MUSCLE) for speed (skips the dominant per-domain cost: a MUSCLE3
            profile alignment against a ~1000-sequence reference set).
        """
        super().__init__()
        if not (0.0 <= threshold <= 1.0):
            raise ValueError("threshold must be between 0 and 1")
        if not (isinstance(keep_top, int) and keep_top >= 1):
            raise ValueError("keep_top must be an integer >= 1")

        self.threshold = threshold
        self.keep_top = keep_top
        self.cache_dir = Path(cache_dir) if cache_dir else DEFAULT_CACHE_DIR
        self.force_retrain = force_retrain
        self.use_muscle_fallback = use_muscle_fallback
        self._model = None

    def _get_model(self):
        if self._model is None:
            self._model = train_model(cache_dir=self.cache_dir, force=self.force_retrain)
        return self._model

    def _results_for(self, domain_id: str, predictions: list) -> list[InferenceResult]:
        unknown = self.result(label="unknown", score=0.0, metadata={})

        results = []
        for pred in predictions:
            if pred.score >= self.threshold:
                metadata = {"smiles": pred.smiles} if pred.smiles else {}
                if not pred.smiles:
                    log.warning(f"no SMILES found for predicted label '{pred.label}'; returning label only")
                results.append(self.result(label=pred.label, score=round(pred.score, 4), metadata=metadata))
            else:
                log.debug(f"prediction '{pred.label}' for sequence {domain_id} below threshold ({pred.score:.4f} < {self.threshold}); skipping")
        return results or [unknown]

    def predict(self, domain: Domain) -> list[InferenceResult]:
        if domain.type != "AMP-binding":
            return []

        model = self._get_model()
        domain_predictions = predict_domains(
            model, {domain.id: domain.sequence}, keep_top=self.keep_top, use_muscle_fallback=self.use_muscle_fallback
        )

        match domain_predictions:
            case []:
                log.warning(f"no A domains found in sequence {domain.id}; unable to predict substrate")
                return [self.result(label="unknown", score=0.0, metadata={})]
            case [a_domain]:
                return self._results_for(domain.id, a_domain.predictions)
            case _:
                log.error(f"found multiple ({len(domain_predictions)}) A domains in sequence {domain.id}; unable to predict substrate")
                return [self.result(label="unknown", score=0.0, metadata={})]

    def predict_many(self, domains: list[Domain]) -> dict[str, list[InferenceResult]]:
        """
        Predict every AMP-binding domain in `domains` in ONE HMMER2/HMMER3(/MUSCLE3
        fallback) pass, instead of one subprocess pair per domain (see
        retromol_antismash.inference.registry.annotate_region, which calls this in
        preference to predict() when it's available -- a region with several NRPS
        modules would otherwise spawn hmmpfam2+hmmscan once *per module*, each
        reloading the full HMM profile database from scratch, which is the dominant
        cost of this CLI-tool-based reimplementation).

        :param domains: candidate domains -- non-"AMP-binding" ones are ignored.
        :return: {domain.id: results}, same per-domain result shape predict() returns.
        """
        amp_domains = [d for d in domains if d.type == "AMP-binding"]
        if not amp_domains:
            return {}

        model = self._get_model()
        sequences = {d.id: d.sequence for d in amp_domains}
        predictions_by_id = {
            p.protein_name: p
            for p in predict_domains(model, sequences, keep_top=self.keep_top, use_muscle_fallback=self.use_muscle_fallback)
        }

        results: dict[str, list[InferenceResult]] = {}
        for d in amp_domains:
            prediction = predictions_by_id.get(d.id)
            if prediction is None:
                log.warning(f"no A domains found in sequence {d.id}; unable to predict substrate")
                results[d.id] = [self.result(label="unknown", score=0.0, metadata={})]
            else:
                results[d.id] = self._results_for(d.id, prediction.predictions)

        return results
