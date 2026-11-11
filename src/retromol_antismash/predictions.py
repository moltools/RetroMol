"""
Resolves substrate/stereochemistry predictions into mxn.yml monomer tokens,
driven entirely by `retromol.data.pmp.yml`.

Two kinds of prediction source are supported, both described generically in
pmp.yml so adding a method never requires touching this module:

- "qualifier": read a value straight out of an antiSMASH GenBank domain's
  raw qualifiers (e.g. `/specificity="substrate consensus: mmal"`).
- "model": read the top-scoring result of a registered `DomainInferenceModel`
  (see retromol_antismash/inference/registry.py) off a domain's annotations
  -- this is how the PARAS NRPS predictor plugs in, and how anyone else's
  custom substrate-prediction model would too.

See pmp.yml's header for the full schema.
"""

from dataclasses import dataclass
from importlib.resources import files
from pathlib import Path
from typing import Any, Literal

import yaml

import retromol.data
from retromol_antismash.model import Domain, InferenceResult


PredictionSource = Literal["qualifier", "model"]


@dataclass(frozen=True)
class PredictionMethod:
    """
    One named prediction method for one axis (e.g. "pks.extender.antismash").

    :param name: the method's name within its axis (the key under `methods.<group>.<axis>` in pmp.yml).
    :param source: "qualifier" to read a GenBank domain qualifier, or "model" to read a registered model's prediction.
    :param combine_with: "reduction_letter" to append the mapped value to a structurally-derived
        PK_A/B/C/D letter, or None if the mapped value is a complete mxn.yml token on its own.
    :param mapping: predicted-value -> mxn.yml token/name. Empty means try the predicted value as a
        literal mxn.yml name, with no rewording.
    :param feature: (qualifier source) antiSMASH feature type the domain lives on, e.g. "aSDomain".
    :param domain: (qualifier source) the domain's `/aSDomain=` value to match, e.g. "PKS_AT".
    :param qualifier: (qualifier source) the GenBank qualifier to read, e.g. "specificity".
    :param key: (qualifier source) the text before ": " in the qualifier line to select.
    :param model: (model source) the registered DomainInferenceModel.name to read predictions from.
    """

    name: str
    source: PredictionSource
    combine_with: str | None
    mapping: dict[str, str]
    feature: str | None = None
    domain: str | None = None
    qualifier: str | None = None
    key: str | None = None
    model: str | None = None

    @classmethod
    def from_dict(cls, name: str, data: dict[str, Any]) -> "PredictionMethod":
        """
        Build a PredictionMethod from its pmp.yml entry.

        :param name: the method's name (the YAML key it was defined under).
        :param data: the method's YAML mapping.
        :return: a PredictionMethod.
        :raises ValueError: if required fields for the method's source are missing.
        """
        source = data.get("source", "qualifier")
        if source == "qualifier":
            missing = [k for k in ("feature", "domain", "qualifier", "key") if k not in data]
            if missing:
                raise ValueError(f"prediction method {name!r} (source=qualifier) missing {missing}")
        elif source == "model":
            if "model" not in data:
                raise ValueError(f"prediction method {name!r} (source=model) missing 'model'")
        else:
            raise ValueError(f"prediction method {name!r} has unknown source {source!r}")

        return cls(
            name=name,
            source=source,
            combine_with=data.get("combine_with"),
            mapping=dict(data.get("mapping") or {}),
            feature=data.get("feature"),
            domain=data.get("domain"),
            qualifier=data.get("qualifier"),
            key=data.get("key"),
            model=data.get("model"),
        )


class PredictionConfig:
    """
    All prediction methods and the active-method selection loaded from pmp.yml.
    """

    def __init__(
        self,
        predictors: dict[str, dict[str, str]],
        methods: dict[str, dict[str, dict[str, PredictionMethod]]],
    ) -> None:
        self.predictors = predictors
        self.methods = methods

    def get_method(self, axis: str) -> PredictionMethod:
        """
        Resolve the active PredictionMethod for an axis.

        :param axis: dotted axis path, e.g. "pks.extender" or "nrps.a_domain".
        :return: the PredictionMethod currently selected for that axis in `predictors`.
        :raises KeyError: if the axis has no selection, or the selection names an undefined method.
        """
        group, _, name = axis.partition(".")
        try:
            selected = self.predictors[group][name]
        except KeyError as exc:
            raise KeyError(f"no predictor selected for {axis!r} in pmp.yml's 'predictors' section") from exc
        try:
            return self.methods[group][name][selected]
        except KeyError as exc:
            raise KeyError(f"predictors.{axis} selects {selected!r}, but methods.{axis}.{selected} isn't defined") from exc

    @classmethod
    def load_from_file(cls, path: str | Path | None = None) -> "PredictionConfig":
        """
        Load a PredictionConfig from a pmp.yml file.

        :param path: path to a pmp.yml file, or None to load the packaged default.
        :return: a PredictionConfig.
        """
        resolved = Path(path) if path else Path(files(retromol.data).joinpath("pmp.yml"))

        with open(resolved, "r") as fo:
            raw = yaml.safe_load(fo) or {}

        predictors = {group: dict(sel or {}) for group, sel in (raw.get("predictors") or {}).items()}

        methods: dict[str, dict[str, dict[str, PredictionMethod]]] = {}
        for group, axes in (raw.get("methods") or {}).items():
            methods[group] = {}
            for axis_name, named_methods in (axes or {}).items():
                methods[group][axis_name] = {
                    method_name: PredictionMethod.from_dict(method_name, method_data)
                    for method_name, method_data in (named_methods or {}).items()
                }

        return cls(predictors, methods)

    @classmethod
    def load_default(cls) -> "PredictionConfig":
        """ Load the packaged default pmp.yml. """
        return cls.load_from_file()


def extract_specificity_value(domain: Domain, qualifier: str, key: str) -> str | None:
    """
    Pull one named value out of a repeated antiSMASH qualifier (e.g.
    `/specificity="substrate consensus: mmal"` -> "mmal" for
    qualifier="specificity", key="substrate consensus").

    :param domain: the Domain whose raw_qualifiers to search.
    :param qualifier: the qualifier name to look in (e.g. "specificity").
    :param key: the text before ": " to match.
    :return: the value after ": ", or None if not present.
    """
    values = domain.raw_qualifiers.get(qualifier)
    if not values:
        return None
    if isinstance(values, str):
        values = [values]

    prefix = f"{key}: "
    for entry in values:
        if entry.startswith(prefix):
            return entry[len(prefix):].strip()
    return None


def resolve_qualifier_method(domain: Domain, method: PredictionMethod) -> str | None:
    """
    Resolve a `source: qualifier` PredictionMethod against a Domain.

    :param domain: the Domain to read the qualifier from.
    :param method: a PredictionMethod with source="qualifier".
    :return: the mapped token fragment, the raw predicted value if `mapping`
        is empty, or None if the qualifier/key wasn't present.
    """
    value = extract_specificity_value(domain, method.qualifier, method.key)
    if value is None:
        return None
    if not method.mapping:
        return value
    return method.mapping.get(value, value)


def resolve_model_result(domain: Domain, method: PredictionMethod) -> InferenceResult | None:
    """
    Get the top-scoring InferenceResult from a `source: model` PredictionMethod's registered model.

    :param domain: the Domain whose annotations to search.
    :param method: a PredictionMethod with source="model".
    :return: the highest-scoring InferenceResult from that model, or None if it produced none.
    """
    if not domain.annotations:
        return None
    results = domain.annotations.by_model(method.model)
    if not results:
        return None
    return max(results, key=lambda r: r.score if r.score is not None else 0.0)


def resolve_model_label(domain: Domain, method: PredictionMethod) -> str | None:
    """
    Resolve a `source: model` PredictionMethod's predicted label, applying
    any name override in `mapping`. Unmapped labels are returned unchanged
    (a caller may then try them as a literal mxn.yml name, or fall back to
    whatever structural resolution the model's own predictions support).

    :param domain: the Domain whose annotations to search.
    :param method: a PredictionMethod with source="model".
    :return: the (possibly remapped) predicted label, or None if the model produced no result.
    """
    result = resolve_model_result(domain, method)
    if result is None:
        return None
    if not method.mapping:
        return result.label
    return method.mapping.get(result.label, result.label)
