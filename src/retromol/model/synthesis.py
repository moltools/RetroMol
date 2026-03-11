"""Data structures for representing synthesis extraction results."""

import logging
from dataclasses import dataclass

from retromol.model.reaction_graph import ReactionGraph


log = logging.getLogger(__name__)


@dataclass(frozen=True)
class SynthesisExtractResult:
    """
    Result of synthesis subgraph extraction.
    
    :var graph: The extracted synthesis subgraph.
    :var solved: Whether the root was successfully solved.
    :var total_cost: Total cost of the extracted subgraph.
    """

    graph: ReactionGraph
    solved: bool
    total_cost: float
