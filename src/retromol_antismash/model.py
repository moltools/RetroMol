"""Module describes the data model for antiSMASH regions and genes."""

from enum import Enum
from dataclasses import dataclass, field, asdict
from typing import Any


class TargetType(Enum):
    """
    Enumeration of possible target types for model inference.
    
    :cvar DOMAIN: Target type representing a domain-level prediction.
    :cvar GENE: Target type representing a gene-level prediction.
    """
    
    DOMAIN = "domain"
    GENE = "gene"


@dataclass
class InferenceResult:
    """
    Represents the result of a model inference.
    
    :param model: The name of the model used for inference.
    :param target: The target type of the inference (e.g., DOMAIN, GENE).
    :param label: The predicted label or class.
    :param score: The confidence score of the prediction (optional).
    :param metadata: Additional metadata related to the inference (optional).
    """

    model: str
    target: TargetType
    label: str
    score: float | None = None
    metadata: dict[str, Any] | None = None

    def to_dict(self) -> dict[str, Any]:
        """
        Convert the InferenceResult instance to a dictionary.
        
        :return: A dictionary representation of the InferenceResult.
        """
        result_dict = asdict(self)
        result_dict["target"] = self.target.value
        return result_dict
    
    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "InferenceResult":
        """
        Create an InferenceResult instance from a dictionary.
        
        :param data: A dictionary containing the inference result data.
        :return: An InferenceResult instance.
        """
        data = data.copy()
        data["target"] = TargetType(data["target"])
        return cls(**data)


@dataclass
class AnnotationSet:
    """
    Represents a set of annotations derived from model inferences.
    
    :param results: A list of inference results.
    """
    results: list[InferenceResult] = field(default_factory=list)

    def add(self, result: InferenceResult) -> None:
        """
        Adds an inference result to the annotation set.
        
        :param result: the inference result to add
        """
        self.results.append(result)

    def by_model(self, model_name: str) -> list[InferenceResult]:
        """
        Retrieves all inference results from a specific model.
        
        :param model_name: The name of the model.
        :return: A list of inference results from the specified model.
        """
        return [r for r in self.results if r.model == model_name]
    
    def to_dict(self) -> dict:
        """
        Converts the annotation set to a dictionary representation.
        
        :return: A dictionary representation of the annotation set.
        """
        return {"results": [r.to_dict() for r in self.results]}
    
    @classmethod
    def from_dict(cls, data: dict) -> "AnnotationSet":
        """
        Creates an AnnotationSet instance from a dictionary representation.
        
        :param data: A dictionary representation of the annotation set.
        :return: An AnnotationSet instance.
        """
        results = [InferenceResult(**r) for r in data.get("results", [])]
        return cls(results=results)


@dataclass
class Domain:
    """
    Represents a biological domain within a gene.
    
    :param id: The unique identifier of the domain.
    :param type: The type or name of the domain.
    :param start: The starting position of the domain within the gene.
    :param end: The ending position of the domain within the gene.
    :param sequence: The amino acid sequence of the domain.
    :param raw_qualifiers: Raw qualifiers or metadata associated with the domain.
    :param annotations: The set of annotations associated with the domain.
    """

    id: str
    type: str
    start: int
    end: int
    sequence: str
    raw_qualifiers: dict[str, Any] = field(default_factory=dict)

    annotations: AnnotationSet = field(default_factory=AnnotationSet)

    def to_dict(self) -> dict:
        """
        Converts the domain to a dictionary representation.
        
        :return: A dictionary representation of the domain.
        """
        data = asdict(self)
        data["annotations"] = self.annotations.to_dict()
        return data
    
    @classmethod
    def from_dict(cls, data: dict) -> "Domain":
        """
        Creates a Domain instance from a dictionary representation.
        
        :param data: A dictionary representation of the domain.
        :return: A Domain instance.
        """
        annotations = AnnotationSet.from_dict(data.get("annotations", {}))
        return cls(
            id=data["id"],
            type=data["type"],
            start=data["start"],
            end=data["end"],
            sequence=data["sequence"],
            raw_qualifiers=data.get("raw_qualifiers", {}),
            annotations=annotations,
        )


class Strand(Enum):
    """
    Represents the strand orientation of a gene.
    
    :cvar FORWARD: The forward strand.
    :cvar REVERSE: The reverse strand.
    """
    FORWARD = "+"
    REVERSE = "-"


@dataclass
class Gene:
    """
    Represents a gene within a biological sequence.
    
    :param id: The unique identifier of the gene.
    :param start: The starting position of the gene within the sequence.
    :param end: The ending position of the gene within the sequence.
    :param strand: The strand orientation of the gene.
    :param sequence: The nucleotide sequence of the gene.
    :param domains: The list of domains associated with the gene.
    :param annotations: the set of annotations associated with the gene
    """

    id: str
    start: int
    end: int
    strand: Strand
    sequence: str

    domains: list[Domain] = field(default_factory=list)
    annotations: AnnotationSet = field(default_factory=AnnotationSet)

    def iter_domains(self) -> list[Domain]:
        """
        Returns the list of domains sorted by their starting position.
        
        :return: A list of domains sorted by start position.
        """
        return sorted(self.domains, key=lambda d: d.start)
    
    def to_dict(self) -> dict:
        """
        Converts the gene to a dictionary representation.
        
        :return: A dictionary representation of the gene.
        """
        data = asdict(self)
        data["strand"] = self.strand.value
        data["domains"] = [d.to_dict() for d in self.domains]
        data["annotations"] = self.annotations.to_dict()
        return data
    
    @classmethod
    def from_dict(cls, data: dict) -> "Gene":
        """
        Creates a Gene instance from a dictionary representation.
        
        :param data: A dictionary representation of the gene.
        :return: A Gene instance.
        """
        domains = [Domain.from_dict(d) for d in data.get("domains", [])]
        annotations = AnnotationSet.from_dict(data.get("annotations", {}))
        return cls(
            id=data["id"],
            start=data["start"],
            end=data["end"],
            strand=Strand(data["strand"]),
            sequence=data["sequence"],
            domains=domains,
            annotations=annotations,
        )


@dataclass
class Region:
    """
    Represents a biological region containing multiple genes.
    
    :param id: The unique identifier of the region.
    :param file_name: The name of the file from which the region was parsed.
    :param start: The starting position of the region within the sequence.
    :param end: The ending position of the region within the sequence.
    :param qualifiers: Additional metadata or qualifiers associated with the region.
    :param genes: The list of genes contained within the region.
    """

    id: str
    file_name: str
    start: int
    end: int
    qualifiers: dict[str, Any] = field(default_factory=dict)
    genes: list[Gene] = field(default_factory=list)

    def iter_genes(self) -> list[Gene]:
        """
        Returns the list of genes sorted by their starting position.
        
        :return: A list of genes sorted by start position.
        """
        return sorted(self.genes, key=lambda g: g.start)

    def to_dict(self) -> dict[str, Any]:
        """
        Converts the Region instance to a dictionary.
        
        :return: A dictionary representation of the Region.
        """
        region_dict = asdict(self)
        region_dict["genes"] = [gene.to_dict() for gene in self.genes]
        return region_dict
    
    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "Region":
        """
        Creates a Region instance from a dictionary.
        
        :param data: A dictionary containing region data.
        :return: A Region instance.
        """
        genes = [Gene.from_dict(gene_data) for gene_data in data.get("genes", [])]
        return cls(
            id=data["id"],
            file_name=data["file_name"],
            start=data["start"],
            end=data["end"],
            qualifiers=data.get("qualifiers", {}),
            genes=genes,
        )
