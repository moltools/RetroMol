#!/usr/bin/env python3
"""
Robust Zenodo uploader:
- Always resolves the *draft deposition* created by /actions/newversion
- Never uploads to the read-only *record* bucket
- Discards stale drafts that block newversion
- Deletes any same-named file in the draft before re-upload
- Publishes after metadata update

Env:
  ZENODO_TOKEN           (required)
  ZENODO_API_URL         (optional; default https://zenodo.org/api)
  ZENODO_DEPOSITION_ID   (required)  # concept/record id you're versioning (e.g., 17410940)
  ZENODO_TITLE           (optional)
  ZENODO_DESCRIPTION     (optional)
  ZENODO_AUTHORS         (optional JSON list, e.g. [{"name":"...","affiliation":"..."}])
  ZENODO_COMMUNITIES     (optional JSON list, e.g. [{"identifier":"..."}])
  ZENODO_KEYWORDS        (optional JSON list of strings)
  VERSION                (required, e.g. v1.2.3 or 1.2.3)
  FILEPATH               (required, path to .zip to upload)
"""

import json
import os
import sys
import pathlib
import requests

API = os.getenv("ZENODO_API_URL", "https://zenodo.org/api").rstrip("/")
TOKEN = os.environ["ZENODO_TOKEN"]
PARAMS = {"access_token": TOKEN}
HEADERS_JSON = {"Content-Type": "application/json"}


def req(method: str, url: str, **kw) -> requests.Response:
    r = requests.request(method, url, **kw)
    if not r.ok:
        # surface server message to logs
        try:
            print(r.text, file=sys.stderr)
        except Exception:
            pass
        r.raise_for_status()
    return r


def get_json(url: str) -> dict:
    return req("GET", url, params=PARAMS).json()


def discard_draft_if_present(concept_id: str) -> None:
    """If a draft exists for this concept, discard it (to avoid 'remove files first')."""
    dep = get_json(f"{API}/deposit/depositions/{concept_id}")
    links = dep.get("links", {})
    discard_url = links.get("discard")
    latest_draft = links.get("latest_draft")
    if latest_draft and discard_url:
        try:
            req("POST", discard_url, params=PARAMS)
        except Exception:
            # If discard not allowed (e.g., already clean), ignore
            pass


def create_new_draft_from(concept_id: str) -> dict:
    """POST newversion and follow Location to the *new draft deposition* JSON."""
    r = requests.post(
        f"{API}/deposit/depositions/{concept_id}/actions/newversion",
        params=PARAMS,
    )
    if r.status_code == 400 and "remove all files" in r.text.lower():
        # Stale draft with files — discard once then retry
        discard_draft_if_present(concept_id)
        r = requests.post(
            f"{API}/deposit/depositions/{concept_id}/actions/newversion",
            params=PARAMS,
        )
    if r.status_code not in (200, 201):
        print(r.text, file=sys.stderr)
        r.raise_for_status()

    loc = r.headers.get("Location")
    if not loc:
        raise RuntimeError("Zenodo did not return a Location header for new draft deposition")
    return get_json(loc)


def ensure_file_not_present_in_draft(draft: dict, filename: str) -> None:
    """If the draft already has a file with this name, delete it so we can re-upload."""
    files = draft.get("files") or []
    for f in files:
        if f.get("filename") == filename:
            file_self = f.get("links", {}).get("self")
            if file_self:
                req("DELETE", file_self, params=PARAMS)


def main() -> None:
    version = os.environ["VERSION"].lstrip("v")  # accept v1.2.3 or 1.2.3
    filepath = pathlib.Path(os.environ["FILEPATH"]).resolve()

    title = os.getenv("ZENODO_TITLE") or f"{pathlib.Path.cwd().name} v{version}"
    description = os.getenv("ZENODO_DESCRIPTION", "")
    authors = json.loads(os.getenv("ZENODO_AUTHORS", "[]"))
    communities = json.loads(os.getenv("ZENODO_COMMUNITIES", "[]"))
    keywords = json.loads(os.getenv("ZENODO_KEYWORDS", "[]"))

    concept_id = os.getenv("ZENODO_DEPOSITION_ID")
    if not concept_id:
        raise SystemExit("ZENODO_DEPOSITION_ID is required (the concept/record id to version).")

    # 1) Make a clean new draft from the concept id and resolve its *draft deposition* JSON
    draft = create_new_draft_from(concept_id)

    # 2) Upload to the *draft deposition's* bucket (never the record bucket)
    links = draft.get("links", {})
    bucket_url = links.get("bucket")
    if not bucket_url:
        # If we somehow didn't get a draft JSON, surface links to debug and abort
        print(json.dumps(links, indent=2), file=sys.stderr)
        raise RuntimeError("Draft deposition has no 'bucket' link; aborting.")

    # Ensure same-named file is not present (idempotent re-runs)
    ensure_file_not_present_in_draft(draft, filepath.name)

    with open(filepath, "rb") as f:
        req("PUT", f"{bucket_url}/{filepath.name}", params=PARAMS, data=f)

    # 3) Update metadata on the draft deposition
    metadata = {
        "metadata": {
            "title": title,
            "upload_type": "software",
            "version": version,
            "description": description or f"Release {version}",
        }
    }
    if authors:
        metadata["metadata"]["creators"] = authors
    if communities:
        metadata["metadata"]["communities"] = communities
    if keywords:
        metadata["metadata"]["keywords"] = keywords

    req("PUT", links["self"], params=PARAMS, data=json.dumps(metadata), headers=HEADERS_JSON)

    # 4) Publish
    req("POST", links["publish"], params=PARAMS)
    print(f"Published Zenodo deposition for v{version} with file {filepath.name}")


if __name__ == "__main__":
    main()
