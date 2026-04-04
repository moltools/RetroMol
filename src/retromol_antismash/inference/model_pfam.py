"""Module for PFAM domain inference model."""

import logging
import math
import os
from typing import TYPE_CHECKING

try:
    from pyhmmer import easel, plan7, hmmer
    PYHMMER_AVAILABLE = True
except ImportError:
    PYHMMER_AVAILABLE = False
    HMM_DB = None

if TYPE_CHECKING:
    from pyhmmer import easel, plan7, hmmer

from retromol_antismash.model import Gene
from retromol_antismash.inference.base import GeneInferenceModel, InferenceResult


log = logging.getLogger(__name__)


def evalue_to_score(evalue: float, cap: float = 50.0) -> float:
    """
    Convert E-value to [0, 1] confidence score.

    :param evalue: E-value from HMMER
    :param cap: value at which score is capped to 1.0
    :return: confidence score between 0 and 1
    """
    if evalue <= 0:
        return 1.0

    score = -math.log10(evalue)

    return min(score / cap, 1.0)


class PfamModel(GeneInferenceModel):
    """
    PFAM domain inference model.

    :param hmm_path: path to the PFAM HMM file
    """

    evalue_cutoff: float = 1e-5
    use_gathering_cutoff: bool = True

    _alphabet: easel.Alphabet | None = None
    _pipeline: plan7.Pipeline | None = None
    _hmm: plan7.HMM | None = None

    def __init__(self, hmm_path: str, label: str) -> None:
        """
        Initialize the PFAM model.

        :param hmm_path: path to the PFAM HMM file
        :param label: label to emit on hit
        """
        super().__init__()

        self.name = f"pfam_{label}"

        self.hmm_path = hmm_path
        self.label = label

    def __post_init__(self) -> None:
        """
        Post-initialization checks.

        :raises FileNotFoundError: if the HMM file does not exist
        """
        super().__post_init__()

        # Make sure hmm_path is valid
        if not os.path.isfile(self.hmm_path):
            log.error(f"PFAM HMM file not found at {self.hmm_path}")
            raise FileNotFoundError(f"PFAM HMM file not found at {self.hmm_path}")
        
    def _init_hmmer(self) -> None:
        """
        Initialize HMMER components.

        :raises ValueError: if no HMM is found in the specified file
        """
        if not PYHMMER_AVAILABLE:
            raise ImportError("HMMER is not available. Please install pyhmmer to use the PFAM model.")

        if self._hmm is not None:
            return
        
        self._alphabet = easel.Alphabet.amino()
        self._pipeline = plan7.Pipeline(self._alphabet, bit_cutoffs="gathering" if self.use_gathering_cutoff else None)

        with plan7.HMMFile(self.hmm_path) as hmm_file:
            self._hmm = list(hmm_file)
        
    def predict(self, gene: Gene) -> list[InferenceResult]:
        """
        Make predictions for a given gene.

        :param gene: the gene to make predictions for
        :return: a list of InferenceResult objects containing the predictions.
        """
        self._init_hmmer()

        protein = gene.sequence
        if not protein:
            return []
        
        text_seq = easel.TextSequence(name=gene.id.encode(), sequence=protein)
        dig_seq = text_seq.digitize(self._alphabet)
        
        hits_iter = hmmer.hmmscan([dig_seq], self._hmm, cpus=1, E=self.evalue_cutoff)

        query_hits = next(hits_iter)  # only one sequence
    
        for hit in query_hits:
            evalue = hit.evalue
            score = evalue_to_score(evalue)
            return [self.result(label=self.label, score=score, metadata={"evalue": evalue})]
        
        return []
