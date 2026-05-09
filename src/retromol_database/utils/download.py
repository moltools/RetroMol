"""Utility functions for downloading and preparing files."""

import os
import shutil
import tarfile
import tempfile
import zipfile
import gzip
from hashlib import sha256
from pathlib import Path
from urllib.request import Request, urlopen

from tqdm import tqdm


def _slug_url(url: str) -> str:
    """
    Deterministic slug from URL (first 16 hex chars of SHA256).

    :param url: URL string
    :return: slug string
    .. note:: this is used to create unique cache directories per URL
    """
    h = sha256(url.encode("utf-8")).hexdigest()[:16]

    # Keep a tiny hint of the basename to help humans
    tail = (_guess_filename_from_url(url) or "payload").split("/")[-1]
    tail = "".join(ch for ch in tail if ch.isalnum() or ch in ("-", "_", "."))

    return f"{tail}-{h}"


def _guess_filename_from_url(url: str) -> str | None:
    """
    Lightweight filename guess from URL path which falls back to None if no sensible name is present.

    :param url: URL string
    :return: filename string or None
    """
    from urllib.parse import urlparse

    path = urlparse(url).path

    if not path:
        return None

    name = Path(path).name

    return name or None


def _download_with_progress(url: str, dest: Path, chunk_size: int = 1024 * 1024) -> None:
    """
    Stream download to `dest` with tqdm progress.

    :param url: URL to download
    :param dest: destination Path
    :param chunk_size: size of read chunks in bytes
    .. note::
        behavior:
        - writes to a temporary .part file and renames on success
        - cleans up partials on error
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp_fd, tmp_path = tempfile.mkstemp(dir=dest.parent, prefix=dest.name + ".", suffix=".part")
    os.close(tmp_fd)

    tmp = Path(tmp_path)
    try:
        req = Request(url, headers={"User-Agent": "biocracker-downloader/1.0"})

        with urlopen(req) as r:
            total = int(r.headers.get("Content-Length") or 0)
            with (
                open(tmp, "wb") as f,
                tqdm(
                    total=total if total > 0 else None,
                    unit="B",
                    unit_scale=True,
                    unit_divisor=1024,
                    desc=f"Downloading {dest.name}",
                ) as pbar,
            ):
                while True:
                    chunk = r.read(chunk_size)
                    if not chunk:
                        break
                    f.write(chunk)
                    pbar.update(len(chunk))

        tmp.replace(dest)

    except Exception:
        # Clean up partials on failure
        try:
            if tmp.exists():
                tmp.unlink()
        finally:
            raise


def _collapse_singleton(path: Path) -> Path:
    """
    Recursively descend through single-child directories.

    :param path: starting Path
    :return: collapsed Path
    .. note::
        stops when:
        - the path is a file, or
        - a directory has 0 or >1 children, i.e., not a singleton
    """
    current = path

    while current.is_dir():
        # Ignore common junk files when deciding "singleton-ness"
        entries = [p for p in current.iterdir() if p.name not in {".DS_Store", "Thumbs.db"}]
        if len(entries) == 1 and entries[0].is_dir():
            current = entries[0]
        else:
            break

    return current


def _resolve_return_path(item_dir: Path) -> Path:
    """
    Choose the most convenient return object:
      - if extracted/ exists:
          * recursively collapse single-directory trees
          * if the final directory contains exactly one file (recursively), return that file
          * else return the (possibly collapsed) directory
      - else (non-archive):
          * if exactly one file in item_dir (ignoring markers), return it
          * if exactly one dir, collapse it using the same rules as above
          * else return item_dir

    :param item_dir: Path to the item directory
    :return: resolved Path
    """
    marker_names = {"URL.txt", ".READY", ".biocracker_cache_marker"}
    extracted = item_dir / "extracted"

    if extracted.exists():
        collapsed = _collapse_singleton(extracted)
        files = [p for p in collapsed.rglob("*") if p.is_file()] if collapsed.is_dir() else [collapsed]
        if len(files) == 1:
            return files[0]
        return collapsed

    # Non-archive case
    items = [p for p in item_dir.iterdir() if p.name not in marker_names]
    files = [p for p in items if p.is_file()]
    dirs = [p for p in items if p.is_dir()]

    if len(files) == 1 and not dirs:
        return files[0]

    if len(dirs) == 1 and not files:
        collapsed = _collapse_singleton(dirs[0])
        files = [p for p in collapsed.rglob("*") if p.is_file()] if collapsed.is_dir() else [collapsed]
        if len(files) == 1:
            return files[0]
        return collapsed

    return item_dir


def download_and_prepare(url: str, base_cache: Path, *, force: bool = False) -> Path:
    """
    Download a URL with a progress bar into the BioCracker cache, delete the temp
    download after completion, unzip if it's a ZIP, and normalize the return:
      - if the result is a single file: return the file path
      - if it's a directory with one file inside: return that file path
      - if multiple files: return the directory path

    :param url: the download URL
    :param base_cache: base cache path
    :param force: if True, re-download and re-prepare even if a ready copy exists
    :return: Path to a file (if resolved to a single file) or a directory
    .. note:: function is idempotent unless `force=True`
    .. note:: handles ZIP archives; other archive formats are not auto-detected
    .. note:: intermediate files are removed on success
    :raises RuntimeError: if extraction fails
    """
    downloads_root = base_cache / "downloads"
    downloads_root.mkdir(parents=True, exist_ok=True)

    slug = _slug_url(url)
    item_dir = downloads_root / slug
    item_dir.mkdir(parents=True, exist_ok=True)

    url_txt = item_dir / "URL.txt"
    ready_marker = item_dir / ".READY"
    archive_path = item_dir / "downloaded"  # temp canonical name; extension added if known

    # Short-circuit if already prepared and not forcing
    if ready_marker.exists() and not force:
        return _resolve_return_path(item_dir)

    # If an archive exists but not ready, try to proceed to preparation
    # Otherwise, (re)download
    if not ready_marker.exists() or force:
        # Decide final download filename (preserve extension if present in URL)
        guessed_name = _guess_filename_from_url(url)
        if guessed_name:
            archive_path = item_dir / guessed_name

        # Only download when missing or force=True
        need_download = force or not archive_path.exists()
        if need_download:
            _download_with_progress(url, archive_path)

        extract_dir = item_dir / "extracted"

        is_zip = zipfile.is_zipfile(archive_path)
        is_tar = tarfile.is_tarfile(archive_path)
        is_gz = archive_path.suffix.lower() == ".gz" and not is_tar

        if is_zip or is_tar or is_gz:
            if extract_dir.exists() and force:
                shutil.rmtree(extract_dir, ignore_errors=True)
            extract_dir.mkdir(parents=True, exist_ok=True)

            if is_gz:
                # gunzip
                out_path = extract_dir / archive_path.with_suffix("").name
                with gzip.open(archive_path, "rb") as fin, open(out_path, "wb") as fout:
                    shutil.copyfileobj(fin, fout, length=1024 * 1024)

            else:
                try:
                    shutil.unpack_archive(str(archive_path), str(extract_dir))
                except shutil.ReadError:
                    if is_zip:
                        with zipfile.ZipFile(archive_path, "r") as zf:
                            zf.extractall(extract_dir)
                    elif is_tar:
                        with tarfile.open(archive_path, "r:*") as tf:
                            tf.extractall(extract_dir)
                    else:
                        raise RuntimeError("unrecognized archive format despite prior checks")

        else:
            # Non-archive: keep the file
            pass

        # Mark URL and READY
        url_txt.write_text(url)
        ready_marker.touch()

    return _resolve_return_path(item_dir)