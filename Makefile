.PHONY: bump-major bump-minor bump-patch build clean release

BUMP?=patch

bump-major:
	python scripts/bump_version.py major
bump-minor:
	python scripts/bump_version.py minor
bump-patch:
	python scripts/bump_version.py patch

build:
	python -m pip install --upgrade build
	python -m build

clean:
	rm -rf dist build .pytest_cache *.egg-info

release:
	gh workflow run release.yml -f bump=$(BUMP)