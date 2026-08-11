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
from urllib.parse import urlsplit

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


def _url_filename(url: str) -> str:
    """Best-effort original filename from a URL's path (ignores query strings)."""
    name = Path(urlsplit(url).path).name
    return name or "download"


def extract(path: str | Path, out_dir: str | Path, source_name: str | None = None) -> None:
    """
    Extract/decompress a downloaded file, or copy it through unchanged.

    :param path: the downloaded file on disk -- often a fixed, extension-less name
        (e.g. "download.raw") the caller chose for the raw download, so it can't be
        used to decide *how* to extract or to name a plain copy/decompressed output.
    :param out_dir: directory to extract/copy into.
    :param source_name: the original filename (e.g. from the source URL), used both
        to pick the extraction strategy and, for the plain-copy/gzip cases, to name
        the resulting file so its extension survives (e.g. "*.sdf" stays findable).
        Falls back to `path`'s own name if not given.
    """
    path = Path(path)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    display_name = source_name or path.name
    name = display_name.lower()

    if name.endswith((".tar.gz", ".tgz", ".tar")):
        mode = "r:gz" if name.endswith((".tar.gz", ".tgz")) else "r"
        with tarfile.open(path, mode) as tf:
            tf.extractall(out_dir)
    elif name.endswith(".zip"):
        with zipfile.ZipFile(path) as zf:
            zf.extractall(out_dir)
    elif name.endswith(".gz"):
        out_path = out_dir / Path(display_name).with_suffix("").name
        with gzip.open(path, "rb") as fin, open(out_path, "wb") as fout:
            shutil.copyfileobj(fin, fout)
    else:
        shutil.copy2(path, out_dir / display_name)


def run(url: str, download_path: str | Path, extract_dir: str | Path | None = None) -> None:
    download(url, download_path)
    if extract_dir is not None:
        extract(download_path, extract_dir, source_name=_url_filename(url))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--url", required=True)
    ap.add_argument("--download-path", required=True, help="where the raw downloaded file is stored")
    ap.add_argument("--extract-dir", default=None, help="if set, extract/decompress the download into this directory")
    args = ap.parse_args()

    run(url=args.url, download_path=args.download_path, extract_dir=args.extract_dir)


if __name__ == "__main__":
    main()
