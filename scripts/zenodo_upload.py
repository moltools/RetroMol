#!/usr/bin/env python3

"""
Upload a file to Zenodo and publish.
Env vars required:
    ZENODO_TOKEN            (required)
    ZENODO_API_URL          (optional; default https://zenodo.org/api/)  # use https://sandbox.zenodo.org/api/ for testing
    ZENODO_DEPOSITION_ID    (optional; if provided, create new version from this deposition)
    ZENODO_TITLE            (optional; defaults to repo name + version)
    ZENODO_DESCRIPTION      (optional)
    ZENODO_AUTHORS          (optional; JSON list of creators: [{"name":"...","affiliation":"..."}])
    ZENODO_COMMUNITIES      (optional; JSON list: [{"identifier":"..."}]) 
    ZENODO_KEYWORDS         (optional; JSON list of strings)
    VERSION                 (required; e.g., 1.2.3)
    FILEPATH                (required; path to file to upload)
"""

import json, os, sys, requests, pathlib

API = os.getenv("ZENODO_API_URL", "https://zenodo.org/api")
TOKEN = os.environ["ZENODO_TOKEN"]
HEADERS = {"Content-Type":"application/json"}
PARAMS = {"access_token": TOKEN}


def req(method: str, url: str, **kw) -> requests.Response:
    """
    Make a request to the Zenodo API and handle errors.
    
    :param method: HTTP method
    :param url: URL to request
    :return: requests.Response object
    """
    r = requests.request(method, url, **kw)
    if not r.ok:
        print(r.text, file=sys.stderr)
        r.raise_for_status()
    return r


def main() -> None:
    """
    Main function to upload and publish file to Zenodo.
    """
    version = os.environ["VERSION"]
    filepath = pathlib.Path(os.environ["FILEPATH"]).resolve()
    title = os.getenv("ZENODO_TITLE") or f"{pathlib.Path.cwd().name} v{version}"
    description = os.getenv("ZENODO_DESCRIPTION", "")
    authors = json.loads(os.getenv("ZENODO_AUTHORS","[]"))
    communities = json.loads(os.getenv("ZENODO_COMMUNITIES","[]"))
    keywords = json.loads(os.getenv("ZENODO_KEYWORDS","[]"))

    dep_id = os.getenv("ZENODO_DEPOSITION_ID")
    if dep_id:
        # Ensure no stuck draft
        existing = req("GET", f"{API}/deposit/depositions/{dep_id}", params=PARAMS).json()
        if "latest_draft" in existing.get("links", {}):
            # discard any existing draft to avoid "Please remove all files first."
            discard_url = existing["links"]["discard"]
            try:
                req("POST", discard_url, params=PARAMS)
                print("Discarded existing draft.")
            except Exception as e:
                print(f"Note: discard failed or not needed: {e}", file=sys.stderr)

        # 2) Create a new version; follow Location to the new deposition
        r = requests.post(f"{API}/deposit/depositions/{dep_id}/actions/newversion", params=PARAMS)
        if r.status_code not in (200, 201):
            print(r.text, file=sys.stderr)
            r.raise_for_status()
        new_draft_url = r.headers.get("Location")
        if not new_draft_url:
            raise RuntimeError("Zenodo did not return Location for new draft.")
        draft = req("GET", new_draft_url, params=PARAMS).json()
    else:
        payload = {"metadata": {"title": title, "upload_type": "software", "version": version}}
        draft = req("POST", f"{API}/deposit/depositions", params=PARAMS,
                    data=json.dumps(payload), headers=HEADERS).json()
    
    print(json.dumps(draft["links"], indent=2))
    bucket_url = draft["links"]["bucket"]
    print("Uploading to bucket:", bucket_url)

    # Upload file into the draft's bucket
    with open(filepath, "rb") as f:
        up = req("PUT", f"{bucket_url}/{filepath.name}", params=PARAMS, data=f)

    # Update metadata (title/version + optional fields)
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

    req("PUT", draft["links"]["self"], params=PARAMS, data=json.dumps(metadata), headers=HEADERS)

    # Publish the deposition
    req("POST", draft["links"]["publish"], params=PARAMS)
    print(f"Published Zenodo deposition for v{version} with file {filepath.name}")


if __name__ == "__main__":
    main()