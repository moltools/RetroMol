"""NPClassifier API client for chemical-class annotation.

Classifies a compound's chemical class directly from its own structure via GNPS2's
free, public NPClassifier service (https://npclassifier.gnps2.org/classify?smiles=...) --
independent of MIBiG's coarse biosynthetic-class labels (PKS/NRPS/...), and available
for NPAtlas compounds too, which MIBiG's biosyn_class never was. A classification has
three list-valued levels (pathway/superclass/class, most to least general) plus a single
is_glycoside boolean; response shape confirmed live against the API (2026-08-26):
{"pathway_results": [...], "superclass_results": [...], "class_results": [...], "isglycoside": bool}.

No published rate limit exists for this service, so callers are responsible for pacing
requests themselves (see annotate_npclassifier.py's --requests-per-second) -- this module
only wraps a single call with retry/backoff.
"""

from __future__ import annotations

import json
import logging
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass

log = logging.getLogger(__name__)

NPCLASSIFIER_URL = "https://npclassifier.gnps2.org/classify"


@dataclass(frozen=True)
class ClassificationResult:
    pathway: list[str]
    superclass: list[str]
    class_: list[str]
    is_glycoside: bool


def classify_smiles(
    smiles: str,
    *,
    timeout: float = 30.0,
    max_retries: int = 3,
    backoff_seconds: float = 2.0,
) -> ClassificationResult | None:
    """Classify one SMILES via NPClassifier, retrying transient failures (timeouts, 5xx,
    429, malformed JSON) with linear backoff. Returns None (logged, not raised) if every
    attempt fails -- one unclassifiable/unreachable-service molecule shouldn't abort a
    whole pipeline run."""
    url = f"{NPCLASSIFIER_URL}?smiles={urllib.parse.quote(smiles, safe='')}"

    for attempt in range(1, max_retries + 1):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as resp:
                data = json.loads(resp.read())
            return ClassificationResult(
                pathway=[str(p) for p in data.get("pathway_results") or []],
                superclass=[str(s) for s in data.get("superclass_results") or []],
                class_=[str(c) for c in data.get("class_results") or []],
                is_glycoside=bool(data.get("isglycoside", False)),
            )
        except (urllib.error.URLError, urllib.error.HTTPError, json.JSONDecodeError, TimeoutError, OSError) as exc:
            if isinstance(exc, urllib.error.HTTPError) and exc.code == 429:
                # Rate-limited: back off much harder than a transient failure, and
                # respect Retry-After if the server sent one -- retrying at the same
                # pace that just got us 429'd only makes the storm worse.
                retry_after = exc.headers.get("Retry-After") if exc.headers else None
                try:
                    delay = float(retry_after) if retry_after is not None else backoff_seconds * attempt * 5
                except ValueError:
                    delay = backoff_seconds * attempt * 5
                log.warning(
                    "NPClassifier rate-limited (429) on attempt %d/%d for %r -- backing off %.1fs",
                    attempt, max_retries, smiles, delay,
                )
            else:
                delay = backoff_seconds * attempt
                log.warning("NPClassifier request failed (attempt %d/%d) for %r: %s", attempt, max_retries, smiles, exc)

            if attempt < max_retries:
                time.sleep(delay)

    log.error("NPClassifier: giving up on %r after %d attempts", smiles, max_retries)
    return None
