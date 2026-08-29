"""
Module for constructing linear readouts from genomic regions.

Note: upstream/downstream scans are genomic (coordinate-based), not biosynthetic!
"""

from abc import ABC, abstractmethod
from collections import Counter
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Literal, overload

from retromol.chem.mol import smiles_to_mol
from retromol.model.rules import RuleSet
from retromol_antismash.model import Region, Gene, Strand, Domain
from retromol_antismash.predictions import (
    PredictionConfig,
    extract_specificity_value,
    resolve_model_result,
    resolve_qualifier_method,
)


DH_TYPES = {"PKS_DH", "PKS_DHt", "PKS_DH2"}
KR_TYPES = {"PKS_KR"}
ER_TYPES = {"PKS_ER"}
KS_TYPES = {"PKS_KS"}
AT_TYPES = {"PKS_AT"}
PKS_TYPES = KS_TYPES | AT_TYPES | KR_TYPES | DH_TYPES | ER_TYPES
PKS_TE_ALIASES = {"Thioesterase", "PKS_TE", "TE"}
PKS_ACCESSORY = KR_TYPES | DH_TYPES | ER_TYPES
PKS_ANCHOR = KS_TYPES


# Common NRPS domain labels found in antiSMASH outputs
NRPS_A = "AMP-binding"
NRPS_C = "Condensation"
NRPS_T_ALIASES = {"PCP", "Thiolation", "T", "Peptidyl-carrier-protein"}
NRPS_E = "Epimerization"
NRPS_MT_ALIASES = {"N-Methyltransferase", "MT"}
NRPS_OX_ALIASES = {"Oxidase", "Ox", "Oxidoreductase"}
NRPS_R_ALIASES = {"Thioester-reductase", "R", "Reductase"}
NRPS_TE = "Thioesterase"


@dataclass(frozen=True)
class DomainRef:
    """
    Reference to a domain within a gene.

    :param gene: Gene object containing the domain.
    :param domain: Domain object within the gene.
    """

    gene: Gene
    domain: Domain


class ModuleType(Enum):
    """
    Enumeration of module types.
    
    :cvar NRPS: Nonribosomal Peptide Synthetase module.
    :cvar PKS: Polyketide Synthase module.
    """

    NRPS = "NRPS"
    PKS = "PKS"


@dataclass
class Module(ABC):
    """
    Base class for a module in a linear readout.

    :param module_index_in_gene: Index of the module within its gene.
    :param start: Starting position of the module.
    :param end: Ending position of the module.
    :param gene_id: ID of the gene containing the module.
    :param gene_strand: Strand of the gene containing the module.
    :param present_domains: List of domain types present in the module.
    """
    module_index_in_gene: int
    start: int
    end: int
    gene_id: str
    gene_strand: Strand
    present_domains: list[str]

    @property
    @abstractmethod
    def type(self) -> ModuleType:
        """
        Abstract property to get the type of the module.
        """
        raise NotImplementedError

    @property
    @abstractmethod
    def substrate(self) -> Any:
        """
        Abstract property to get the substrate information for the module.
        """
        raise NotImplementedError
    
    def to_dict(self) -> dict[str, Any]:
        """
        Convert the Module object to a dictionary representation.

        :return: Dictionary representation of the Module.
        """
        raise NotImplementedError
    
    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "Module":
        """
        Create a Module object from a dictionary representation.

        :param data: Dictionary representation of the Module.
        :return: Module object.
        """
        raise NotImplementedError
    

@dataclass
class NRPSAnatomy:
    """
    Anatomy of a Nonribosomal Peptide Synthetase (NRPS) module.

    :param has_C: Presence of condensation domain.
    :param has_T: Presence of thiolation domain.
    :param has_E: Presence of epimerization domain.
    :param has_MT: Presence of methyltransferase domain.
    :param has_Ox: Presence of oxidase domain.
    :param has_R: Presence of reductase domain.
    :param has_TE: Presence of thioesterase domain.
    """

    has_C: bool
    has_T: bool
    has_E: bool
    has_MT: bool
    has_Ox: bool
    has_R: bool
    has_TE: bool

    def to_dict(self) -> dict[str, bool]:
        """
        Convert the NRPSAnatomy object to a dictionary representation.

        :return: Dictionary representation of the NRPSAnatomy.
        """
        return {
            "has_C": self.has_C,
            "has_T": self.has_T,
            "has_E": self.has_E,
            "has_MT": self.has_MT,
            "has_Ox": self.has_Ox,
            "has_R": self.has_R,
            "has_TE": self.has_TE,
        }
    
    @classmethod
    def from_dict(cls, data: dict[str, bool]) -> "NRPSAnatomy":
        """
        Create a NRPSAnatomy object from a dictionary representation.

        :param data: Dictionary representation of the NRPSAnatomy.
        :return: NRPSAnatomy object.
        """
        return cls(
            has_C=data.get("has_C", False),
            has_T=data.get("has_T", False),
            has_E=data.get("has_E", False),
            has_MT=data.get("has_MT", False),
            has_Ox=data.get("has_Ox", False),
            has_R=data.get("has_R", False),
            has_TE=data.get("has_TE", False),
        )


@dataclass
class NRPSSubstrate:
    """
    Substrate information for a Nonribosomal Peptide Synthetase (NRPS) module.

    :param name: Name of the predicted substrate.
    :param smiles: SMILES representation of the substrate.
    :param score: Confidence score of the substrate prediction.
    """

    name: str | None
    smiles: str | None
    score: float | None

    def to_dict(self) -> dict[str, Any]:
        """
        Convert the NRPSSubstrate object to a dictionary representation.

        :return: Dictionary representation of the NRPSSubstrate.
        """
        return {
            "name": self.name,
            "smiles": self.smiles,
            "score": self.score,
        }
    
    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "NRPSSubstrate":
        """
        Create a NRPSSubstrate object from a dictionary representation.

        :param data: Dictionary representation of the NRPSSubstrate.
        :return: NRPSSubstrate object.
        """
        return cls(
            name=data.get("name", None),
            smiles=data.get("smiles", None),
            score=data.get("score", None),
        )


class ATLoadingMode(Enum):
    """
    Enumeration of acyltransferase (AT) loading modes.

    :cvar CIS: Cis-acting AT domain.
    :cvar TRANS: Trans-acting AT domain.
    :cvar UNKNOWN: Unknown AT loading mode.
    """

    CIS = "cis"
    TRANS = "trans"
    UNKNOWN = "unknown"


@dataclass 
class PKSAnatomy:
    """
    Anatomy of a Polyketide Synthase (PKS) module.

    :param has_active_KR: Presence of active ketoreductase domain.
    :param has_active_DH: Presence of active dehydratase domain.
    :param has_active_ER: Presence of active enoylreductase domain.
    :param has_AT: Presence of acyltransferase domain.
    :param extender_digit: Extender-unit substituent digit resolved from the AT domain's
        substrate-specificity prediction (see pmp.yml's "pks.extender" axis), e.g. "2" for
        methylmalonyl-CoA. None when unresolved (e.g. trans-AT modules, or an unmapped
        antiSMASH substrate call).
    :param beta_stereo: R/S suffix for the beta-hydroxyl carbon, resolved from the KR domain's
        stereochemistry prediction (see pmp.yml's "pks.kr_stereochemistry" axis, "stereo_suffix"
        combine mode). None when unresolved or (by design) when that axis's mapping is empty.
    """
    AT_loading_mode: ATLoadingMode

    has_active_KR: bool
    has_active_DH: bool
    has_active_ER: bool
    extender_digit: str | None = None
    beta_stereo: str | None = None

    def to_dict(self) -> dict[str, Any]:
        """
        Convert the PKSAnatomy object to a dictionary representation.

        :return: Dictionary representation of the PKSAnatomy.
        """
        return {
            "AT_loading_mode": self.AT_loading_mode.value,
            "has_active_KR": self.has_active_KR,
            "has_active_DH": self.has_active_DH,
            "has_active_ER": self.has_active_ER,
            "extender_digit": self.extender_digit,
            "beta_stereo": self.beta_stereo,
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "PKSAnatomy":
        """
        Create a PKSAnatomy object from a dictionary representation.

        :param data: Dictionary representation of the PKSAnatomy.
        :return: PKSAnatomy object.
        """
        raw = (data.get("AT_loading_mode") or "unknown").lower()
        AT_loading_mode = ATLoadingMode(raw) if raw in {"cis","trans","unknown"} else ATLoadingMode.UNKNOWN

        return cls(
            AT_loading_mode=AT_loading_mode,
            has_active_KR=data.get("has_active_KR", False),
            has_active_DH=data.get("has_active_DH", False),
            has_active_ER=data.get("has_active_ER", False),
            extender_digit=data.get("extender_digit", None),
            beta_stereo=data.get("beta_stereo", None),
        )


class PKSExtenderUnit(Enum):
    """
    Enumeration of PKS extender unit types.

    :cvar PKS_A: PKS extender unit type A.
    :cvar PKS_B: PKS extender unit type B.
    :cvar PKS_C: PKS extender unit type C.
    :cvar PKS_D: PKS extender unit type D.
    :cvar UNCLASSIFIED: Unclassified extender unit type.
    """

    PKS_A = "PK_A"
    PKS_B = "PK_B"
    PKS_C = "PK_C"
    PKS_D = "PK_D"
    UNCLASSIFIED = "UNCLASSIFIED"


@dataclass
class PKSSubstrate:
    """
    Substrate information for a Polyketide Synthase (PKS) module.

    :param extender_unit: Type of extender unit used in the PKS module.
    :param substituent_type: Extender-unit substituent digit (e.g. "2" for methylmalonyl-derived),
        resolved via pmp.yml's "pks.extender" prediction axis. None when unresolved.
    :param beta_stereo: R/S suffix for the beta-hydroxyl carbon, resolved via pmp.yml's
        "pks.kr_stereochemistry" prediction axis. None when unresolved.
    """

    extender_unit: PKSExtenderUnit
    substituent_type: str | None = None
    beta_stereo: str | None = None

    def to_dict(self) -> dict[str, str]:
        """
        Convert the PKSSubstrate object to a dictionary representation.

        :return: Dictionary representation of the PKSSubstrate.
        """
        return {
            "extender_unit": self.extender_unit.value,
            "substituent_type": self.substituent_type,
            "beta_stereo": self.beta_stereo,
        }

    @classmethod
    def from_dict(cls, data: dict[str, str]) -> "PKSSubstrate":
        """
        Create a PKSSubstrate object from a dictionary representation.

        :param data: Dictionary representation of the PKSSubstrate.
        :return: PKSSubstrate object.
        """
        return cls(
            extender_unit=PKSExtenderUnit(data.get("extender_unit", "UNCLASSIFIED")),
            substituent_type=data.get("substituent_type", None),
            beta_stereo=data.get("beta_stereo", None),
        )


@dataclass
class NRPSModule(Module):
    """
    Nonribosomal peptide synthetase (NRPS) module.

    :param anatomy: Anatomical features of the NRPS module.
    :param substrate: Predicted substrate information for the NRPS module.
    """

    anatomy: NRPSAnatomy    
    predicted_substrate: NRPSSubstrate | None = None

    @property
    def type(self) -> ModuleType:
        """
        Get the type of the module.

        :return: ModuleType.NRPS
        """
        return ModuleType.NRPS

    @property
    def substrate(self) -> NRPSSubstrate | None:
        """
        Get the predicted substrate information for the NRPS module.

        :return: NRPSSubstrate object containing substrate information, or None if not available.
        """
        return self.predicted_substrate
    
    def to_dict(self) -> dict[str, Any]:
        """
        Convert the NRPSModule object to a dictionary representation.

        :return: Dictionary representation of the NRPSModule.
        """
        return {
            "type": self.type.value,
            "module_index_in_gene": self.module_index_in_gene,
            "start": self.start,
            "end": self.end,
            "gene_id": self.gene_id,
            "gene_strand": self.gene_strand.value,
            "present_domains": self.present_domains,
            "anatomy": self.anatomy.to_dict(),
            "predicted_substrate": self.predicted_substrate.to_dict() if self.predicted_substrate else None,
        }
    
    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "NRPSModule":
        """
        Create a NRPSModule object from a dictionary representation.

        :param data: Dictionary representation of the NRPSModule
        :return: NRPSModule object
        """
        anatomy_data = data.get("anatomy", {})
        substrate_data = data.get("predicted_substrate", None)

        return cls(
            module_index_in_gene=data["module_index_in_gene"],
            start=data["start"],
            end=data["end"],
            gene_id=data["gene_id"],
            gene_strand=Strand(data["gene_strand"]),
            present_domains=data["present_domains"],
            anatomy=NRPSAnatomy.from_dict(anatomy_data),
            predicted_substrate=NRPSSubstrate.from_dict(substrate_data) if substrate_data else None,
        )


@dataclass
class PKSModule(Module):
    """
    Polyketide synthase (PKS) module.

    :param type: Module type (PKS)
    :param anatomy: Anatomical features of the PKS module.
    """

    anatomy: PKSAnatomy

    @property
    def type(self) -> ModuleType:
        """
        Get the type of the module.

        :return: ModuleType.PKS
        """
        return ModuleType.PKS

    @property
    def substrate(self) -> PKSSubstrate:
        """
        Get the predicted substrate information for the PKS module.

        :return: PKSSubstrate object containing substrate information.
        """
        # Configure factory type
        def setup_substrate(extender_unit: PKSExtenderUnit) -> PKSSubstrate:
            return PKSSubstrate(
                extender_unit=extender_unit,
                substituent_type=self.anatomy.extender_digit,
                beta_stereo=self.anatomy.beta_stereo,
            )

        # Rules:
        # - KS + AT with neither KR nor DH nor ER => PKS_A
        # - KS + AT + KR (no DH and no ER) => PKS_B (KR after AT is naturally true in window order)
        # - KS + AT + KR + DH (no ER) => PKS_C
        # - KS + AT + KR + DH + ER => PKS_D
        # - else UNCLASSIFIED
        # Note: assumes that presence of AT domain is already established

        # True activity from qualifiers
        KR = self.anatomy.has_active_KR
        DH = self.anatomy.has_active_DH
        ER = self.anatomy.has_active_ER

        # Product state logic
        eff_DH = DH and KR          # DH needs KR product to act in canonical cycle
        eff_ER = ER and KR and DH   # ER typically needs DH product (enoyl)

        match (KR, eff_DH, eff_ER):
            case (False, _,     _    ): return setup_substrate(PKSExtenderUnit.PKS_A)
            case (True,  False, _    ): return setup_substrate(PKSExtenderUnit.PKS_B)
            case (True,  True,  False): return setup_substrate(PKSExtenderUnit.PKS_C)
            case (True,  True,  True ): return setup_substrate(PKSExtenderUnit.PKS_D)
            case _:                     return setup_substrate(PKSExtenderUnit.UNCLASSIFIED)

    def to_dict(self) -> dict[str, Any]:
        """
        Convert the PKSModule object to a dictionary representation.

        :return: Dictionary representation of the PKSModule.
        """
        return {
            "type": self.type.value,
            "module_index_in_gene": self.module_index_in_gene,
            "start": self.start,
            "end": self.end,
            "gene_id": self.gene_id,
            "gene_strand": self.gene_strand.value,
            "present_domains": self.present_domains,
            "anatomy": self.anatomy.to_dict(),
        }
    
    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "PKSModule":
        """
        Create a PKSModule object from a dictionary representation.

        :param data: Dictionary representation of the PKSModule.
        :return: PKSModule object.
        """
        anatomy_data = data.get("anatomy", {})

        return cls(
            module_index_in_gene=data["module_index_in_gene"],
            start=data["start"],
            end=data["end"],
            gene_id=data["gene_id"],
            gene_strand=Strand(data["gene_strand"]),
            present_domains=data["present_domains"],
            anatomy=PKSAnatomy.from_dict(anatomy_data),
        )


@dataclass
class LinearReadout:
    """
    A linear readout consisting of a sequence of modules.
    
    :param id: unique Identifier for the linear readout.
    :param file_name: Source file name from which the linear readout was derived.
    :param start: Starting position of the linear readout.
    :param end: Ending position of the linear readout.
    :param qualifiers: Additional metadata or qualifiers associated with the linear readout.
    :param modules: List of modules in the linear readout.
    """
    
    id: str
    file_name: str
    start: int
    end: int
    qualifiers: dict[str, Any] = field(default_factory=dict)

    modules: list[Module] = field(default_factory=list)
    modifiers: list[str] = field(default_factory=list)

    def __str__(self) -> str:
        """
        String representation of the LinearReadout.
        
        :return: String representation of the LinearReadout.
        """
        return f"LinearReadout(id={self.id}, start={self.start}, end={self.end}, modules={len(self.modules)})"
    
    @property
    def num_modules(self) -> int:
        """
        Get the number of modules in the linear readout.

        :return: Number of modules.
        """
        return len(self.modules)
    
    def to_dict(self) -> dict[str, Any]:
        """
        Convert the LinearReadout object to a dictionary representation.

        :return: Dictionary representation of the LinearReadout.
        """
        return {
            "id": self.id,
            "file_name": self.file_name,
            "start": self.start,
            "end": self.end,
            "qualifiers": self.qualifiers,
            "modules": [module.to_dict() for module in self.modules],
            "modifiers": self.modifiers,
        }
    
    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "LinearReadout":
        """
        Create a LinearReadout object from a dictionary representation.

        :param data: Dictionary representation of the LinearReadout.
        :return: LinearReadout object.
        """
        modules_data = data.get("modules", [])
        modules: list[Module] = []

        for mod_data in modules_data:
            mod_type = mod_data.get("type", None)
            if mod_type == ModuleType.NRPS.value:
                modules.append(NRPSModule.from_dict(mod_data))
            elif mod_type == ModuleType.PKS.value:
                modules.append(PKSModule.from_dict(mod_data))
            else:
                raise ValueError(f"Unknown module type: {mod_type}")

        return cls(
            id=data["id"],
            file_name=data["file_name"],
            start=data["start"],
            end=data["end"],
            qualifiers=data.get("qualifiers", {}),
            modules=modules,
            modifiers=data.get("modifiers", []),
        )
    
    @overload
    def biosynthetic_order(self, by_orf: Literal[False] = False) -> list[Module]: ...
    @overload
    def biosynthetic_order(self, by_orf: Literal[True] = True) -> list[tuple[str, list[Module]]]: ...

    def biosynthetic_order(self, by_orf: bool = False):
        """
        Return modules in biosynthetic order.

        :param by_orf: If True, group modules by their originating gene (ORF).
        :return: List of Module objects in biosynthetic order, or list of tuples (gene_id, list of Module) if by_orf is True.
        """
        if not self.modules:
            return []
        
        # Group modules by gene
        by_gene: dict[str, list[Module]] = {}
        for m in self.modules:
            by_gene.setdefault(m.gene_id, []).append(m)

        # Infer strand per gene (sanity check)
        gene_strand: dict[str, Strand] = {}
        for gid, mods in by_gene.items():
            s = mods[0].gene_strand
            if any(m.gene_strand is not s for m in mods):
                raise ValueError(f"mixed gene_strand in gene_id={gid}")
            gene_strand[gid] = s

        # Infer global biosyntehtic direction
        strand_counts = Counter(gene_strand.values())
        global_reverse = strand_counts[Strand.REVERSE] > strand_counts[Strand.FORWARD]

        # Order genes along biosynthetic direction using genomic position
        gene_ids = sorted(
            by_gene.keys(),
            key=lambda gid: min(m.start for m in by_gene[gid]),
            reverse=global_reverse,
        )

        if by_orf:
            grouped: list[tuple[str, list[Module]]] = []
            for gid in gene_ids:
                mods = by_gene[gid]
                if gene_strand[gid] is Strand.FORWARD:
                    mods_sorted = sorted(mods, key=lambda m: m.start)
                else:
                    mods_sorted = sorted(mods, key=lambda m: m.start, reverse=True)
                grouped.append((gid, mods_sorted))

            return grouped

        # Flatten modules in biosynthetic order
        out: list[Module] = []
        for gid in gene_ids:
            mods = by_gene[gid]
            if gene_strand[gid] is Strand.FORWARD:
                mods_sorted = sorted(mods, key=lambda m: m.start)
            else:
                mods_sorted = sorted(mods, key=lambda m: m.start, reverse=True)
            out.extend(mods_sorted)

        return out


def _domain_types(domains: list[Domain]) -> set[str]:
    """
    Helper function to extract the set of domain types from a list of Domain objects.
    
    :param domains: List of Domain objects.
    :return: Set of domain type strings.
    """
    return {d.type for d in domains if d.type is not None}


def _is_domain_type(domain: Domain, label: str | set[str]) -> bool:
    """
    Check if a domain matches a given type label or set of labels.

    :param domain: Domain object to check.
    :param label: Domain type label or set of labels to match against.
    :return: True if the domain type matches the label(s), False otherwise.
    """
    if not domain.type:
        return False

    if isinstance(label, set):
        return domain.type in label
    
    return domain.type == label


def _is_pks_ks(d: Domain) -> bool:
    """
    Check if a domain is a PKS KS domain.
    
    :param d: Domain object to check.
    :return: True if the domain is a PKS KS domain, False otherwise.
    """
    return d.type == "PKS_KS"


def _is_pks_domain(d: Domain) -> bool:
    """
    Check if a domain is a PKS domain.
    
    :param d: Domain object to check.
    :return: True if the domain is a PKS domain, False otherwise.
    """
    return d.type in PKS_TYPES or (d.type in PKS_TE_ALIASES)
 

def _is_active_accessory_domain(domain: Domain) -> bool:
    """
    Determine if an accessory domain (KR, DH, ER) is active based on its qualifiers.
    
    :param domain: Domain object to evaluate.
    :return: True if the domain is active, False if inactive.
    """
    if not domain.type:
        return True  # can't tell, assume active
    
    if domain.type not in PKS_ACCESSORY:
        return True  # not a reducible domain, consider active by default
    
    texts = []
    if domain.id:
        texts.append(domain.id)
    for _, vals in domain.raw_qualifiers.items():
        if isinstance(vals, (list, tuple)):
            texts.extend(map(str, vals))
        else:
            texts.append(str(vals))

    blob = " ".join(texts).lower()

    # Common antiSMASH phrasing patterns
    inactive_flags = [
        "inactive",
        "nonfunctional",
        "non-functional",
        "inactivated",
        "broken",
        "truncated",
    ]
    is_active = not any(flag in blob for flag in inactive_flags)

    return is_active


def _classify_pks_window(window: list[Domain]) -> tuple[set[str], bool, bool, bool, bool]:
    """
    Classify a PKS module window based on the presence and activity of domains.

    :param window: List of Domain objects in the module window.
    :return: tuple containing:
        - module type (str)
        - set of present domain types (set[str])
        - has active KR (bool)
        - has active DH (bool)
        - has active ER (bool)
        - has AT (bool)
    """
    types_linear = [d.type for d in window if d.type in PKS_TYPES]
    present = set(types_linear)

    has_AT = "PKS_AT" in present
    has_active_KR = any(d.type in KR_TYPES and _is_active_accessory_domain(d) for d in window)
    has_active_DH = any(d.type in DH_TYPES and _is_active_accessory_domain(d) for d in window)
    has_active_ER = any(d.type in ER_TYPES and _is_active_accessory_domain(d) for d in window)

    return present, has_active_KR, has_active_DH, has_active_ER, has_AT


def _is_AT_only_gene(gene: Gene) -> bool:
    """
    Helper function to determine if a gene is an acyltransferase-domain-only gene.
    
    :param g: Gene object.
    :return: True if the gene is an AT-only gene, False otherwise.
    """
    types = _domain_types(gene.domains) 
    return ("PKS_AT" in types) and all(t in {"PKS_AT"} for t in types)


def _find_genomic_upstream_AT_only_gene(all_genes: list[Gene], gene_idx_in_genomic_order: int) -> Gene | None:
    """
    Return the nearest upstream gene that is AT-only (relative to all_genes order).

    :param all_genes: List of Gene objects.
    :param gene_idx_in_genomic_order: Index of the current gene in all_genes.
    :return: Gene object of the nearest upstream AT-only gene, or None if not found.
    """
    for j in range(gene_idx_in_genomic_order - 1, -1, -1):
        if _is_AT_only_gene(all_genes[j]):
            return all_genes[j]
        
    return None


def genes_biosynthetic(region: Region) -> list[Gene]:
    """
    Return genes in biosynthetic order within a region.

    :param region: Region object.
    :return: List of Gene objects in biosynthetic order.
    """
    genes = list(region.iter_genes())
    strand_counts = Counter(g.strand for g in genes)
    global_reverse = strand_counts[Strand.REVERSE] > strand_counts[Strand.FORWARD]
    return sorted(genes, key=lambda g: g.start, reverse=global_reverse)


def domains_biosynthetic(gene: Gene) -> list[Domain]:
    """
    Return domains in biosynthetic order within a gene.
    
    :param gene: Gene object.
    :return: List of Domain objects in biosynthetic order.
    .. note:: We assume Domain.start/end are genomic coordinates.
    """
    doms = sorted(gene.domains, key=lambda d: d.start)
    if gene.strand is Strand.REVERSE:
        doms = list(reversed(doms))

    return doms


def region_domain_stream(region: Region) -> list[DomainRef]:
    """
    Return domains in biosynthetic order within a region.
    
    :param region: Region object.
    :return: List of DomainRef objects in biosynthetic order.
    """
    out: list[DomainRef] = []
    for g in genes_biosynthetic(region):
        for d in domains_biosynthetic(g):
            out.append(DomainRef(gene=g, domain=d))
    
    return out


def collect_nrps_modules(gene: Gene, config: PredictionConfig | None = None) -> list[NRPSModule]:
    """
    Collect NRPS modules from a given gene.

    :param gene: Gene object to analyze.
    :param config: prediction configuration (see pmp.yml) to resolve A-domain substrates with,
        or None to load the packaged default (see pmp.yml's "predictors.nrps.a_domain").
    :return: List of NRPSModule objects.
    """
    config = config or PredictionConfig.load_default()
    a_domain_method = config.get_method("nrps.a_domain")

    doms: list[Domain] = domains_biosynthetic(gene)
    out: list[NRPSModule] = []

    # Indices of A domains in left-to-right order
    a_idx = [i for i, d in enumerate(doms) if _is_domain_type(d, NRPS_A)]
    if not a_idx:
        return out  # no A domains, no modules
    
    for mi, ai in enumerate(a_idx):
        # Extend window backward by one if there is an immediately previous C (same gene)
        start_i = ai
        if ai - 1 >= 0 and _is_domain_type(doms[ai - 1], NRPS_C):
            start_i = ai - 1
        
        # Extend forward until (but not including) the next A-domain
        end_i = a_idx[mi + 1] if mi + 1 < len(a_idx) else len(doms)

        window = doms[start_i:end_i]
        present = _domain_types(window)

        has_C = any(_is_domain_type(d, NRPS_C) for d in window)
        has_T = any(_is_domain_type(d, NRPS_T_ALIASES) for d in window)
        has_E = any(_is_domain_type(d, NRPS_E) for d in window)
        has_MT = any(_is_domain_type(d, NRPS_MT_ALIASES) for d in window)
        has_Ox = any(_is_domain_type(d, NRPS_OX_ALIASES) for d in window)
        has_R = any(_is_domain_type(d, NRPS_R_ALIASES) for d in window)
        has_TE = any(_is_domain_type(d, NRPS_TE) for d in window)

        s = min(d.start for d in window)
        e = max(d.end for d in window)

        # Retrieve A domain substrate specificity prediction, via whichever method
        # pmp.yml's "predictors.nrps.a_domain" currently selects.
        A = doms[ai]
        substrate_pred: NRPSSubstrate | None = None

        if a_domain_method.source == "model":
            top_pred = resolve_model_result(A, a_domain_method)
            if top_pred:
                name = top_pred.label
                if a_domain_method.mapping:
                    name = a_domain_method.mapping.get(name, name)
                substrate_pred = NRPSSubstrate(
                    name=name,
                    smiles=(top_pred.metadata or {}).get("smiles", None),
                    score=top_pred.score,
                )
        elif a_domain_method.source == "qualifier":
            resolved = resolve_qualifier_method(A, a_domain_method)
            if resolved is not None:
                substrate_pred = NRPSSubstrate(name=resolved, smiles=None, score=None)

        out.append(NRPSModule(
            module_index_in_gene=mi,
            start=s,
            end=e,
            gene_id=gene.id,
            gene_strand=gene.strand,
            present_domains=list(present),
            anatomy=NRPSAnatomy(
                has_C=has_C,
                has_T=has_T,
                has_E=has_E,
                has_MT=has_MT,
                has_Ox=has_Ox,
                has_R=has_R,
                has_TE=has_TE,
            ),
            predicted_substrate=substrate_pred,
        ))

    return out


def collect_pks_modules(
    region: Region,
    max_cross_gene_bp: int = 20_000,
    config: PredictionConfig | None = None,
) -> list[PKSModule]:
    """
    Collect PKS modules across a genomic region, allowing for cross-gene module assembly.

    :param region: Region object representing the genomic region.
    :param max_cross_gene_bp: Maximum base pair distance to search across genes for module assembly.
    :param config: prediction configuration (see pmp.yml) to resolve extender-unit substituents
        with, or None to load the packaged default.
    :return: List of PKSModule objects collected across the region.
    """
    config = config or PredictionConfig.load_default()
    extender_method = config.get_method("pks.extender")
    kr_stereo_method = config.get_method("pks.kr_stereochemistry")

    stream = region_domain_stream(region)

    # Locate all KS anchors in the stream
    ks_pos = [i for i, ref in enumerate(stream) if _is_pks_ks(ref.domain)]
    if not ks_pos:
        return []  # no KS domains, no modules
    
    out: list[PKSModule] = []
    module_index_by_gene: dict[str, int] = Counter()

    for k_i, start_idx in enumerate(ks_pos):
        end_idx = ks_pos[k_i + 1] if k_i + 1 < len(ks_pos) else len(stream)
        ks_ref = stream[start_idx]
        ks = ks_ref.domain

        # Cancidate window: KS -> next KS (exclusive)
        window_refs = stream[start_idx:end_idx]

        # Don't vaccum up far-away stuff
        filtered: list[DomainRef] = []
        ks_end = ks.end
        for ref in window_refs:
            d = ref.domain
            if d is ks:
                filtered.append(ref)
                continue
            if abs(d.start - ks_end) <= max_cross_gene_bp:
                filtered.append(ref)
            else:
                # Too far away; stop early
                break
        
        # Collect PKS domains in the window
        window_domains = [r.domain for r in filtered if _is_pks_domain(r.domain)]
        (
            present,
            has_active_KR,
            has_active_DH,
            has_active_ER,
            has_AT
        ) = _classify_pks_window(window_domains)

        # Determine AT mode (cis or trans)
        # Note that upstream should be upstream in gene list here (genomic order)
        genes = list(region.iter_genes())
        gene_idx = genes.index(ks_ref.gene)
        if has_AT:
            AT_src: ATLoadingMode = ATLoadingMode.CIS
        else:
            AT_src = (
                ATLoadingMode.TRANS
                if _find_genomic_upstream_AT_only_gene(genes, gene_idx)
                else ATLoadingMode.UNKNOWN
            )

        # DHt is more commonly found in trans PKS modules, so we treat it as inactive in cis modules
        extender_digit: str | None = None
        if AT_src is ATLoadingMode.CIS:
            present_DH_types = present.intersection(DH_TYPES)
            if len(present_DH_types) == 1 and "PKS_DHt" in present_DH_types:
                has_active_DH = False

            # Resolve the extender-unit substituent digit from the module's own AT domain.
            # Trans-AT modules borrow their AT from elsewhere in the record (or lack one
            # entirely), so this is only attempted for cis-AT modules with a domain in hand.
            at_domain = next((d for d in window_domains if d.type in AT_TYPES), None)
            if at_domain is not None and extender_method.source == "qualifier":
                extender_digit = resolve_qualifier_method(at_domain, extender_method)

        # Resolve the beta-hydroxyl stereo suffix from the module's own KR domain.
        # "stereo_suffix" methods deliberately have no-passthrough semantics: an empty
        # (or non-matching) mapping means "don't apply a suffix", since antiSMASH's raw
        # "A1"/"B2" class labels aren't R/S descriptors on their own (see pmp.yml).
        beta_stereo: str | None = None
        if kr_stereo_method.source == "qualifier" and kr_stereo_method.combine_with == "stereo_suffix":
            kr_domain = next((d for d in window_domains if d.type in KR_TYPES), None)
            if kr_domain is not None:
                raw_kr_class = extract_specificity_value(kr_domain, kr_stereo_method.qualifier, kr_stereo_method.key)
                if raw_kr_class is not None:
                    beta_stereo = kr_stereo_method.mapping.get(raw_kr_class)

        # Use window_domains bounds for start/end
        s = min(r.domain.start for r in filtered)
        e = max(r.domain.end for r in filtered)

        gid = ks_ref.gene.id
        mi = module_index_by_gene[gid]
        module_index_by_gene[gid] += 1

        out.append(PKSModule(
            module_index_in_gene=mi,
            start=s,
            end=e,
            gene_id=gid,
            gene_strand=ks_ref.gene.strand,
            present_domains=list(present),
            anatomy=PKSAnatomy(
                AT_loading_mode=AT_src,
                has_active_KR=has_active_KR,
                has_active_DH=has_active_DH,
                has_active_ER=has_active_ER,
                extender_digit=extender_digit,
                beta_stereo=beta_stereo,
            ),
        ))

    return out


def linear_readout(region: Region, config: PredictionConfig | None = None) -> LinearReadout:
    """
    Construct a linear readout from the given genomic region.

    :param region: Region object representing the genomic region.
    :param config: prediction configuration (see pmp.yml) to resolve PKS/NRPS substrates with,
        or None to load the packaged default.
    :return: LinearReadout object containing the collected modules.
    """
    assert isinstance(region, Region), "region must be an instance of Region"

    config = config or PredictionConfig.load_default()

    collected: list[Module] = []
    modifiers: list[str] = []

    # Collect NRPS modules (gene-level)
    for gene in region.iter_genes():
        collected.extend(collect_nrps_modules(gene, config=config))

    # Collect PKS modules region-wide (cross-gene)
    collected.extend(collect_pks_modules(region, config=config))

    # Check if there are any gene-level modifiers
    for gene in region.iter_genes():
        if gene.annotations:
            for result in gene.annotations.results:
                label = result.label
                modifiers.append(label)

    return LinearReadout(
        id=region.id,
        file_name=region.file_name,
        start=region.start,
        end=region.end,
        qualifiers=region.qualifiers,
        modules=collected,
        modifiers=modifiers,
    )


def module_primary_sequence_tokens(module: Module, ruleset: RuleSet) -> tuple[str, list[str]]:
    """
    Map a single module's predicted substrate onto the ruleset's monomer vocabulary,
    so BGC readouts land in the same name/fingerprint-token space as compound primary
    sequences (see `retromol_synthesis.reconstruction.reconstruct_linear_readout`,
    whose `primary_sequence` entries are matching-rule names drawn from the same
    ruleset).

    The two module types resolve very differently, because they carry different kinds
    of evidence:

    - A PKS module's extender unit is known at the reduction-level for certain (its
      `PKSSubstrate.extender_unit`, e.g. "PK_A"), and additionally at the specific
      side-chain level (`PKSSubstrate.substituent_type`, e.g. "2") whenever pmp.yml's
      "pks.extender" predictor resolved one from the AT domain's substrate-specificity
      call. When a specific name (e.g. "A2") exists as a rule in `ruleset`, it's used;
      otherwise this falls back to the group-level pseudonym tokens every rule at that
      reduction level shares (e.g. ["PK_A", "PK"]) -- naming a specific rule without one
      existing would claim identity information we don't actually have.
    - An NRPS module's `NRPSSubstrate.name` is tried first as a literal
      rule name (this is where pmp.yml's "nrps.a_domain" `mapping` overrides take effect,
      and where a `source: qualifier` fallback with no SMILES gets resolved) -- but before
      that literal lookup, two anatomy-driven name transforms are tried first, most
      specific first, since both are downstream modifications no A-domain substrate
      caller (Stachelhaus code, NRPSPredictor2, PARAS) can see: they all predict identity
      purely from the A-domain's binding-pocket residues, which only decide which side
      chain gets activated, so `substrate.name` is always the plain, unmethylated
      L-residue's identity, regardless of what the rest of the module does to it.
        1. `NRPSAnatomy.has_MT` (an N-methyltransferase domain, whether standalone or
           embedded in the A-domain, present in the module) prefixes the name with
           "N-methyl" (e.g. "leucine" -> "N-methylleucine") -- the MT domain methylates
           the backbone amide nitrogen of the loaded aminoacyl-PCP intermediate, not the
           alpha carbon, so this is independent of stereochemistry.
        2. `NRPSAnatomy.has_E` (an epimerization domain present) appends a "^L"/"^D"
           stereo suffix (D if present, L otherwise) -- see point 1's independence: this
           applies on top of, not instead of, the "N-methyl" prefix.
      Both transforms are tried together first (most specific: "N-methyl<name>^L" or
      "N-methyl<name>^D"), then MT-only ("N-methyl<name>"), then E-only ("<name>^L" or
      "<name>^D"), then the bare `name`, so any level of an unmapped combination degrades
      gracefully instead of returning "X" outright -- most exotic/non-proteinogenic
      monomers, and achiral residues like glycine, don't have every variant defined in
      mxn.yml. Neither transform is a pmp.yml-configurable predictor axis like
      "pks.kr_stereochemistry" -- MT/E domain presence is read directly off the module's
      anatomy, not a qualifier/model prediction, and the "N-methyl"/"^D" mxn.yml token
      construction is a fixed biological fact rather than something a predictor swap
      should override. If nothing in the candidate list matches, and a SMILES is
      available (the default PARAS path), it falls back to matching that SMILES against
      every rule by structure (`RuleSet.find_structural_matches`) -- PARAS' raw labels
      aren't guaranteed to match a rule's `name` string verbatim. Whether that structural
      match also requires matching stereochemistry follows `ruleset.match_stereochemistry`,
      the same toggle every other structural match in this codebase respects -- PARAS'
      predicted substrate is a specific named compound (its SMILES comes from a fixed
      label->structure lookup, not a per-domain stereochemistry guess), so requiring
      an exact stereo match here is meaningful, not just permissively ignored.

    :param module: a PKS or NRPS module from a BGC's LinearReadout.
    :param ruleset: the rule set to resolve the module's substrate against.
    :return: a tuple of (display name, fingerprint tokens) for this module. The
        display name is "X" and tokens are empty whenever nothing could be resolved --
        the same "unidentified block" convention used for compound primary sequences.
    """
    def _by_name(name: str) -> tuple[str, list[str]] | None:
        matches = [r for r in ruleset.matching_rules if r.name == name]
        if not matches:
            return None
        tokens: set[str] = set()
        for rule in matches:
            tokens.add(rule.name)
            tokens.update(rule.pseudonyms)
        return min(rule.name for rule in matches), sorted(tokens)

    if module.type is ModuleType.PKS:
        substrate = module.substrate
        extender_unit = substrate.extender_unit
        if extender_unit is PKSExtenderUnit.UNCLASSIFIED:
            return "X", []

        if substrate.substituent_type:
            letter = extender_unit.value.rsplit("_", 1)[-1]  # "PK_A" -> "A"
            base_name = f"{letter}{substrate.substituent_type}"

            # mxn.yml's beta-hydroxyl stereo suffix sits right after the letter, before
            # the digit (e.g. "B^R2", not "B2^R") -- only meaningful for "B" (the only
            # reduction level with a free beta-hydroxyl: A never has an active KR, C's
            # is consumed by dehydration into an E/Z alkene instead of an R/S center, and
            # D's is consumed by the enoylreductase). This intentionally only reaches the
            # beta-only fallback form (e.g. "B^R2"), not the doubly-stereo-specified forms
            # that also fix the alpha carbon (e.g. "B^R2^R") -- resolving the alpha center
            # isn't wired up yet (see pmp.yml's "pks.kr_stereochemistry" comments).
            if substrate.beta_stereo and letter == "B":
                stereo_specific = _by_name(f"B^{substrate.beta_stereo}{substrate.substituent_type}")
                if stereo_specific is not None:
                    return stereo_specific

            specific = _by_name(base_name)
            if specific is not None:
                return specific

        return extender_unit.value, [extender_unit.value, "PK"]

    # NRPS: try the resolved substrate name with both anatomy-driven transforms applied
    # (an "N-methyl" prefix if `has_MT`, a "^L"/"^D" stereo suffix picked by `has_E`),
    # then each transform alone, then the bare resolved name as a literal rule name
    # (this last one is where an explicit pmp.yml mapping, or a qualifier-sourced
    # antiSMASH call with no SMILES, gets resolved), then fall back to structural SMILES
    # matching. See docstring above for why both transforms have to be tried independent
    # of whether mxn.yml actually has a rule for every combination.
    substrate = module.substrate
    if substrate and substrate.name:
        stereo_suffix = "D" if module.anatomy.has_E else "L"
        methylated_name = f"N-methyl{substrate.name}" if module.anatomy.has_MT else None

        candidates = []
        if methylated_name:
            candidates.append(f"{methylated_name}^{stereo_suffix}")
            candidates.append(methylated_name)
        candidates.append(f"{substrate.name}^{stereo_suffix}")
        candidates.append(substrate.name)

        for candidate in candidates:
            resolved = _by_name(candidate)
            if resolved is not None:
                return resolved

    smiles = substrate.smiles if substrate else None
    if not smiles:
        return "X", []

    try:
        mol = smiles_to_mol(smiles)
    except ValueError:
        return "X", []

    matches = ruleset.find_structural_matches(mol, match_stereochemistry=ruleset.match_stereochemistry)
    if not matches:
        return "X", []

    tokens: set[str] = set()
    for rule in matches:
        tokens.add(rule.name)
        tokens.update(rule.pseudonyms)

    # Deterministic choice among structurally-tied rules (e.g. two rules that
    # describe the same free-molecule structure but differ in which atom continues
    # the polymer chain, such as "aspartic acid" vs. an isoAsp-linked variant).
    name = min(rule.name for rule in matches)

    return name, sorted(tokens)


def bgc_primary_sequence(readout: LinearReadout, ruleset: RuleSet) -> tuple[list[str], list[list[str]]]:
    """
    Convert a BGC's linear readout into a primary sequence representation comparable
    to a compound's: one display name and one fingerprint token list per module, both
    drawn from the same matching-rule vocabulary, in biosynthetic order.

    :param readout: the BGC linear readout to convert.
    :param ruleset: the rule set to resolve each module's substrate against.
    :return: a tuple of (display names, per-module fingerprint tokens), one entry per module in biosynthetic order.
    """
    names: list[str] = []
    tokens: list[list[str]] = []

    for module in readout.biosynthetic_order():
        name, module_tokens = module_primary_sequence_tokens(module, ruleset)
        names.append(name)
        tokens.append(module_tokens)

    return names, tokens