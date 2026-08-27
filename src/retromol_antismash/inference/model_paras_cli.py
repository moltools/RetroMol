"""
Alternative to model_paras.ParasModel: predicts NRPS A-domain substrate specificity using
retromol_paras, a self-contained reimplementation of PARAS that runs domain extraction via
command-line HMMER2/HMMER3/MUSCLE3 (no pyhmmer) and trains its own RandomForestClassifier
against this repo's pinned scikit-learn version (no unpickling someone else's model file).

Registered under a different DomainInferenceModel.name ("paras_cli") than the original
("paras") so both can coexist -- point pmp.yml's "predictors.nrps.a_domain" at "paras_cli"
to use this one (see retromol_paras/__init__.py for the required system binaries).
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
    ) -> None:
        """
        :param threshold: minimum predicted probability for a substrate call to be kept.
        :param keep_top: number of top-scoring substrate predictions to keep per domain.
        :param cache_dir: directory to cache the trained model and extracted training
            signatures in (see retromol_paras.train.train_model). Defaults to ~/.retromol_paras.
        :param force_retrain: retrain from scratch even if a cached model is present.
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
        self._model = None

    def _get_model(self):
        if self._model is None:
            self._model = train_model(cache_dir=self.cache_dir, force=self.force_retrain)
        return self._model

    def predict(self, domain: Domain) -> list[InferenceResult]:
        if domain.type != "AMP-binding":
            return []

        unknown = self.result(label="unknown", score=0.0, metadata={})
        model = self._get_model()

        domain_predictions = predict_domains(
            model, {domain.id: domain.sequence}, keep_top=self.keep_top
        )

        match domain_predictions:
            case []:
                log.warning(f"no A domains found in sequence {domain.id}; unable to predict substrate")
                return [unknown]
            case [a_domain]:
                results = []
                for pred in a_domain.predictions:
                    if pred.score >= self.threshold:
                        metadata = {"smiles": pred.smiles} if pred.smiles else {}
                        if not pred.smiles:
                            log.warning(f"no SMILES found for predicted label '{pred.label}'; returning label only")
                        results.append(self.result(label=pred.label, score=round(pred.score, 4), metadata=metadata))
                    else:
                        log.debug(f"prediction '{pred.label}' for sequence {domain.id} below threshold ({pred.score:.4f} < {self.threshold}); skipping")
                return results or [unknown]
            case _:
                log.error(f"found multiple ({len(domain_predictions)}) A domains in sequence {domain.id}; unable to predict substrate")
                return [unknown]
