"""ChEMBL bulk-database client for bioactivity annotation.

Downloads and extracts ChEMBL's official SQLite release (ftp.ebi.ac.uk) once, then
looks up a compound's known medicinal-relevance signals by its standard InChIKey --
ChEMBL's own maximum clinical-development phase (MAX_PHASE) and WHO ATC therapeutic
category (a defined vocabulary, unlike MIBiG's free-text chem_acts) -- rather than
querying ChEMBL's live API per compound, which would be far too slow at our scale.

Schema confirmed live against ChEMBL's own schema documentation
(https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/schema_documentation.txt,
2026-08-26): COMPOUND_STRUCTURES.STANDARD_INCHI_KEY -> MOLECULE_DICTIONARY.MOLREGNO
(MAX_PHASE, CHEMBL_ID) -> MOLECULE_ATC_CLASSIFICATION -> ATC_CLASSIFICATION.LEVEL1_DESCRIPTION.
"""

from __future__ import annotations

import logging
import sqlite3
import tarfile
import urllib.request
from dataclasses import dataclass
from pathlib import Path

from tqdm import tqdm

log = logging.getLogger(__name__)

# "latest" always resolves to the newest ChEMBL release -- pin to a specific
# release's URL (e.g. .../chembl_37/chembl_37_sqlite.tar.gz) if a pipeline run needs
# to be reproducible across ChEMBL releases.
CHEMBL_SQLITE_URL = "https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/chembl_37_sqlite.tar.gz"

# ChEMBL's MAX_PHASE values (MOLECULE_DICTIONARY.MAX_PHASE) mapped to a display label.
# Only phases with real clinical-development meaning are surfaced -- 0/-1 (no reported
# development / unknown) aren't stored, the same presence-only convention every other
# annotation table uses. Python's dict lookup treats 3 and 3.0 as the same key, so this
# works whether sqlite hands back MAX_PHASE as an int or a float.
MAX_PHASE_LABELS = {
    4: "Approved",
    3: "Phase 3",
    2: "Phase 2",
    1: "Phase 1",
    0.5: "Early phase 1",
}


def download_chembl_sqlite(dest_dir: str | Path, *, url: str = CHEMBL_SQLITE_URL, force: bool = False) -> Path:
    """Download and extract ChEMBL's SQLite release into `dest_dir`, returning the path
    to the extracted .db file (no-op if one is already there, unless `force`).

    This is a multi-GB download (~5.4GB compressed as of ChEMBL 37) -- streamed to disk
    in chunks rather than read into memory, unlike taxonomy.py's much smaller taxdump.
    """
    dest_dir = Path(dest_dir).expanduser()
    dest_dir.mkdir(parents=True, exist_ok=True)

    existing = sorted(dest_dir.rglob("*.db"))
    if not force and existing:
        return existing[0]

    archive_path = dest_dir / "chembl_sqlite.tar.gz"
    log.info("downloading ChEMBL SQLite release from %s", url)
    with urllib.request.urlopen(url) as resp, open(archive_path, "wb") as out:
        total = int(resp.headers.get("Content-Length") or 0) or None
        with tqdm(total=total, desc="download_chembl", unit="B", unit_scale=True, unit_divisor=1024) as pbar:
            while chunk := resp.read(1024 * 1024):
                out.write(chunk)
                pbar.update(len(chunk))

    log.info("extracting %s", archive_path)
    with tarfile.open(archive_path, mode="r:gz") as tar:
        members = tar.getmembers()
        for member in tqdm(members, desc="extract_chembl", unit="file"):
            tar.extract(member, path=dest_dir)
    archive_path.unlink()

    found = sorted(dest_dir.rglob("*.db"))
    if not found:
        raise FileNotFoundError(f"ChEMBL archive at {url} did not contain a .db file")
    return found[0]


@dataclass(frozen=True)
class AtcCategory:
    level1_description: str  # display label -- coarse WHO therapeutic category
    level5: str  # e.g. "J01FA01" -- the compound's own specific ATC code, for linking out


@dataclass(frozen=True)
class ChemblBioactivity:
    chembl_id: str
    max_phase_label: str | None
    atc_categories: list[AtcCategory]


class ChemblDB:
    """Read-only lookup over a local ChEMBL SQLite release, by standard InChIKey."""

    def __init__(self, sqlite_path: str | Path) -> None:
        self.con = sqlite3.connect(f"file:{Path(sqlite_path).expanduser()}?mode=ro", uri=True)

    def close(self) -> None:
        self.con.close()

    def __enter__(self) -> "ChemblDB":
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self.close()

    def bioactivity_for_inchikey(self, inchikey: str) -> ChemblBioactivity | None:
        """None if `inchikey` has no match in ChEMBL at all; otherwise a result whose
        `max_phase_label`/`atc_categories` may still both be empty (a compound can be
        in ChEMBL -- have activity/assay records -- without being an ATC-classified drug
        or having reached any clinical phase)."""
        row = self.con.execute(
            """
            SELECT md.molregno, md.chembl_id, md.max_phase
            FROM compound_structures cs
            JOIN molecule_dictionary md ON md.molregno = cs.molregno
            WHERE cs.standard_inchi_key = ?
            """,
            (inchikey,),
        ).fetchone()
        if row is None:
            return None

        molregno, chembl_id, max_phase = row
        max_phase_label = MAX_PHASE_LABELS.get(max_phase)

        atc_rows = self.con.execute(
            """
            SELECT DISTINCT ac.level1_description, ac.level5
            FROM molecule_atc_classification mac
            JOIN atc_classification ac ON ac.level5 = mac.level5
            WHERE mac.molregno = ?
            """,
            (molregno,),
        ).fetchall()
        atc_categories = [
            AtcCategory(level1_description=str(r[0]), level5=str(r[1])) for r in atc_rows if r[0] and r[1]
        ]

        return ChemblBioactivity(
            chembl_id=str(chembl_id),
            max_phase_label=max_phase_label,
            atc_categories=atc_categories,
        )
