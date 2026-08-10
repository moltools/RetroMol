"""Step 2: download a source file/archive and optionally extract it.

Handles the two shapes the pipeline needs: a single (possibly gzipped) file, like
the NPAtlas SDF, and a tar.gz/zip archive containing many files, like MIBiG's JSON
and GenBank bundles.
"""

import argparse
import gzip
import shutil
import tarfile
import zipfile
from pathlib import Path

import requests


def download(url: str, dest: str | Path) -> None:
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    with requests.get(url, stream=True, timeout=300) as resp:
        resp.raise_for_status()
        with open(dest, "wb") as fh:
            for chunk in resp.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    fh.write(chunk)


def extract(path: str | Path, out_dir: str | Path) -> None:
    path = Path(path)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    name = path.name.lower()
    if name.endswith((".tar.gz", ".tgz", ".tar")):
        mode = "r:gz" if name.endswith((".tar.gz", ".tgz")) else "r"
        with tarfile.open(path, mode) as tf:
            tf.extractall(out_dir)
    elif name.endswith(".zip"):
        with zipfile.ZipFile(path) as zf:
            zf.extractall(out_dir)
    elif name.endswith(".gz"):
        out_path = out_dir / path.with_suffix("").name
        with gzip.open(path, "rb") as fin, open(out_path, "wb") as fout:
            shutil.copyfileobj(fin, fout)
    else:
        shutil.copy2(path, out_dir / path.name)


def run(url: str, download_path: str | Path, extract_dir: str | Path | None = None) -> None:
    download(url, download_path)
    if extract_dir is not None:
        extract(download_path, extract_dir)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--url", required=True)
    ap.add_argument("--download-path", required=True, help="where the raw downloaded file is stored")
    ap.add_argument("--extract-dir", default=None, help="if set, extract/decompress the download into this directory")
    args = ap.parse_args()

    run(url=args.url, download_path=args.download_path, extract_dir=args.extract_dir)


if __name__ == "__main__":
    main()
