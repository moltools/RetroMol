#!/usr/bin/env python3

"""
Bump semantic version in pyproject.toml.

Usage: 
    python scripts/bump_version.py [major|minor|patch]
"""

from os import major
import argparse, re, sys, pathlib

ROOT = pathlib.Path(__file__).resolve().parents[1]
PYPROJECT = ROOT / "pyproject.toml"


def read_pyprojects() -> tuple[str, tuple[int, int, int]]:
    """
    Read pyproject.toml and extract the current version.
    
    :return: tuple of file content and version tuple (major, minor, patch)
    """
    txt = PYPROJECT.read_text(encoding="utf-8")
    m = re.search(r'(?m)^\s*version\s*=\s*"(\d+)\.(\d+)\.(\d+)"\s*$', txt)
    if not m:
        raise SystemExit("Could not find version = \"X.Y.Z\" in pyproject.toml")
    return txt, tuple(map(int, m.groups()))


def write_pyproject(txt: str, new_version: str) -> None:
    """
    Write the new version back to pyproject.toml.

    :param txt: original file content
    :param new_version: new version string
    :return: None
    """
    new_txt = re.sub(
        r'(?m)^(\s*version\s*=\s*")(\d+\.\d+\.\d+)(")',
        r"\g<1>" + new_version + r"\3",
        txt,
        count=1,
    )
    PYPROJECT.write_text(new_txt, encoding="utf-8")


def bump(ver: tuple[int, int, int], kind: str) -> str:
    """
    Bump the version based on the kind.

    :param ver: current version tuple (major, minor, patch)
    :param kind: kind of bump ("major", "minor", "patch")
    :return: new version string
    """
    major, minor, patch = ver
    if kind == "major":
        return f"{major + 1}.0.0"
    elif kind == "minor":
        return f"{major}.{minor + 1}.0"
    elif kind == "patch":
        return f"{major}.{minor}.{patch + 1}"
    else:
        raise SystemExit(f"Kind must be one of major, minor, patch; got {kind!r}")
    

def main() -> None:
    """
    Main function to parse arguments and bump version.
    """
    p = argparse.ArgumentParser()
    p.add_argument("kind", choices=["major", "minor", "patch"])
    args = p.parse_args()

    txt, cur = read_pyprojects()
    new_version = bump(cur, args.kind)
    write_pyproject(txt, new_version)
    print(new_version)  # workflow reads new version from stdout for tagging


if __name__ == "__main__":
    main()
