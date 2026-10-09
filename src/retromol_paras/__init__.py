"""
retromol_paras -- a self-contained reimplementation of PARAS (Adenylation-domain
substrate specificity prediction) for RetroMol.

This does NOT import or depend on the `parasect` package
(https://github.com/BTheDragonMaster/parasect). It vendors just enough of that
project's published logic (domain extraction via HMMER2/HMMER3/MUSCLE3 command-line
tools, extended-signature featurisation, and a scikit-learn RandomForestClassifier
trained on the same public dataset) to reproduce PARAS predictions, while pinning
scikit-learn to a version this repo controls rather than whatever a pretrained
`.paras` model file happened to be pickled with.

Requires the following external command-line tools on PATH (none of them are
pip-installable): `hmmpfam2` (HMMER2), `hmmscan` (HMMER3), `hmmpress` (HMMER3,
used once to prepare the packaged HMM3 profile), and `muscle` (MUSCLE v3 --
the `-in1/-in2/-profile` CLI, NOT MUSCLE v5's incompatible interface).

See environment.yml (repo root) for a conda environment that installs all of
these, including an Apple Silicon fallback if the arm64 builds don't solve or
run.

See retromol_paras.train for building/caching a model, and retromol_paras.predict
for running predictions against it.
"""
