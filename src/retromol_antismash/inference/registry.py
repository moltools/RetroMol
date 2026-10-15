"""Registry for inference models."""

import logging
from dataclasses import dataclass

from retromol_antismash.model import Region
from retromol_antismash.inference.base import GeneInferenceModel, DomainInferenceModel


GENE_MODELS: list[GeneInferenceModel] = []
DOMAIN_MODELS: list[DomainInferenceModel] = []


log = logging.getLogger(__name__)


def register_gene_model(model: GeneInferenceModel) -> None:
    """
    Register a gene inference model.
    
    :param model: the gene inference model to register
    :raises TypeError: if the model is not an instance of GeneInferenceModel
    """
    # Make sure model is an instance of GeneInferenceModel
    if not isinstance(model, GeneInferenceModel):
        raise TypeError("model must be an instance of GeneInferenceModel")

    # Check that the model is not already registered; if it is, do not register it again
    registered = [m for m in GENE_MODELS if m.name == model.name]
    if registered:
        return

    GENE_MODELS.append(model)


def register_domain_model(model: DomainInferenceModel) -> None:
    """
    Register a domain inference model.

    :param model: the domain inference model to register
    """
    # Make sure model is an instance of DomainInferenceModel
    if not isinstance(model, DomainInferenceModel):
        raise TypeError("model must be an instance of DomainInferenceModel")

    # Check that the model is not already registered; if it is, do not register it again
    registered = [m for m in DOMAIN_MODELS if m.name == model.name]
    if registered:
        return

    DOMAIN_MODELS.append(model)


def get_gene_models() -> list[GeneInferenceModel]:
    """
    Get the list of registered gene inference models.
    
    :return: list of registered gene inference models
    """
    return GENE_MODELS


def get_domain_models() -> list[DomainInferenceModel]:
    """
    Get the list of registered domain inference models.

    :return: list of registered domain inference models
    """
    return DOMAIN_MODELS


@dataclass(frozen=True)
class Ctx:
    """
    Context information for logging.
    
    :cvar region: optional region identifier
    :cvar gene: optional gene identifier
    :cvar domain: optional domain identifier
    :cvar model: optional model identifier

    .. note::
        Usage example:
            ctx = Ctx(region="chr1", gene="BRCA1")
            logger.info(f"{ctx.prefix()}This is a log message.")
    """

    region: str | None = None
    gene: str | None = None
    domain: str | None = None
    model: str | None = None

    def prefix(self) -> str:
        """
        Generate a log prefix string based on the context.
        
        :return: formatted prefix string
        """
        parts = []
        if self.region: parts.append(f"region={self.region}")
        if self.gene: parts.append(f"gene={self.gene}")
        if self.domain: parts.append(f"domain={self.domain}")
        if self.model: parts.append(f"model={self.model}")
        return ("[" + " ".join(parts) + "] ") if parts else ""


def annotate_region(
    region: Region,
    gene_models: list[GeneInferenceModel] | None = None,
    domain_models: list[DomainInferenceModel] | None = None,
) -> None:
    """
    Annotate all domains in all genes of the given region.

    :param region: the genomic region to annotate
    :param gene_models: gene inference models to run, or None to use every globally
        registered gene model (`get_gene_models()`).
    :param domain_models: domain inference models to run, or None to use every
        globally registered domain model (`get_domain_models()`). Pass this
        explicitly rather than registering globally when a model needs per-call
        configuration (e.g. PARAS' probability threshold) -- `register_domain_model`
        keeps whichever instance was registered first, so it can't express "use a
        different threshold for this call."
    """
    gene_models = get_gene_models() if gene_models is None else gene_models
    domain_models = get_domain_models() if domain_models is None else domain_models

    log.debug(Ctx(region=region.id).prefix() + f"annotating {len(region.genes)} genes")

    for gene in region.iter_genes():
        gctx = Ctx(region=region.id, gene=gene.id)
        log.debug(gctx.prefix() + "annotating gene")

        # Gene inference
        for m in gene_models:
            mctx = Ctx(region=region.id, gene=gene.id, model=m.name)
            log.debug(mctx.prefix() + "running gene inference")

            results = m.predict(gene)
            for r in results:
                gene.annotations.add(r)

            log.debug(mctx.prefix() + f"added {len(results)} results")

    # Domain inference. Batched per model across the WHOLE region (not per gene, and
    # not per domain) wherever a model supports it -- a model like paras_cli that
    # shells out to command-line tools pays a large fixed cost per call (reloading an
    # entire HMM profile database from scratch), so calling predict() once per domain
    # instead of predict_many() once for every domain a model cares about is a major
    # (not just marginal) slowdown, not a style choice. Models without predict_many
    # (e.g. gene-level-adjacent PFAM domain models) fall back to the original
    # per-domain predict() loop, unaffected.
    all_domains = [domain for gene in region.iter_genes() for domain in gene.iter_domains()]

    for m in domain_models:
        predict_many = getattr(m, "predict_many", None)
        if predict_many is None:
            continue

        mctx = Ctx(region=region.id, model=m.name)
        log.debug(mctx.prefix() + f"running batched domain inference over {len(all_domains)} domains")

        results_by_id = predict_many(all_domains)
        for domain in all_domains:
            for r in results_by_id.get(domain.id, []):
                domain.annotations.add(r)

    unbatched_models = [m for m in domain_models if getattr(m, "predict_many", None) is None]
    if unbatched_models:
        for domain in all_domains:
            dctx = Ctx(region=region.id, domain=domain.id)

            for m in unbatched_models:
                mctx = Ctx(region=region.id, domain=domain.id, model=m.name)
                log.debug(mctx.prefix() + f"running domain inference ({domain.type})")

                results = m.predict(domain)
                for r in results:
                    domain.annotations.add(r)

                log.debug(mctx.prefix() + f"added {len(results)} results")