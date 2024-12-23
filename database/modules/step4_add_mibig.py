import typing as ty
import os
import re
import json
import joblib
from collections import defaultdict
from Bio.Seq import Seq  # biopython
from parasect.api import run_parasect # pip install git+https://github.com/BTheDragonMaster/parasect.git@webapp


TEMP_DIR = os.path.join(os.path.dirname(__file__), "temp")
MODEL_DIR = os.path.join(os.path.dirname(__file__), "models")
PARASECT_MODEL = joblib.load(os.path.join(MODEL_DIR, "model.parasect"))
PARASECT_MODEL.set_params(n_jobs=1)


class Gene:
    def __init__(self, name: str, start: int, end: int, strand: str) -> None:
        self.name = name
        self.start = start
        self.end = end
        self.strand = strand
        self.modules = list()

    def __str__(self) -> str:
        return f"Gene(name={self.name}, start={self.start}, end={self.end}, strand={self.strand})"

    def __hash__(self) -> int:
        return hash(self.name)
    
    def __eq__(self, other: ty.Any) -> bool:
        if not isinstance(other, Gene): return False
        return self.name == other.name
    

class ProtoCluster:
    def __init__(self, record_id: str, start: int, end: int, product: str, category: str) -> None:
        self.record_id = record_id
        self.start = start
        self.end = end
        self.product = product
        self.category = category
        self.genes = list()

    def __str__(self) -> str:
        return f"ProtoCluster(start={self.start}, end={self.end}, product={self.product}, category={self.category})"

    def __hash__(self) -> int:
        return hash((self.start, self.end, self.product, self.category))
    
    def __eq__(self, other: ty.Any) -> bool:
        if not isinstance(other, ProtoCluster): return False
        return (self.start == other.start and self.end == other.end and self.product == other.product and self.category == other.category)
    
    def __len__(self) -> int:
        return len(self.genes)
    
    def to_motif_code(self, dna: str, predict_specificities: bool = True) -> ty.Dict[str, ty.Any]:
        """Convert the proto cluster to a motif code."""
        to_predict = []  # collect sequences for PARAS A-domain specificity prediction

        # bundle relevant modules into putative domains
        groups = []
        for gene in self.genes:
            gene_sequence_added = False
            modules = gene.modules
            current_group = []
            for module in modules:
                if module in ["PKS_KS", "AMP-binding"]:
                    if module == "AMP-binding" and not gene_sequence_added:
                        gene_start = gene.start
                        gene_end = gene.end
                        gene_strand = gene.strand
                        gene_dna = dna[gene_start:gene_end]
                        if gene_strand == "-": gene_dna = str(Seq(gene_dna).reverse_complement())
                        gene_protein = str(Seq(gene_dna).translate())
                        gene_fasta = f">{gene.name}\n{gene_protein}"
                        to_predict.append(gene_fasta)
                        gene_sequence_added = True
                    if current_group: groups.append(current_group)
                    current_group = [module]
                elif (module in ["PKS_DH", "PKS_ER", "PKS_KR"] and current_group and current_group[0] == "PKS_KS"): current_group.append(module)
                else: continue 

            if current_group: groups.append(current_group)

        # predict specificity for A-domains
        to_predict = "\n".join(to_predict)
        if not predict_specificities or not to_predict: predicted_specificities = []
        else:
            try:
                results = run_parasect(selected_input=to_predict, selected_input_type="fasta", path_temp_dir=TEMP_DIR, model=PARASECT_MODEL)
                predicted_specificities = []
                for result in results:
                    result_json = result.to_json()
                    predictions = result_json["predictions"]
                    predictions = sorted(predictions, key=lambda x: x["probability"], reverse=True)
                    top_prediction = predictions[0]
                    top_prediction_name = top_prediction["substrate_name"]
                    predicted_specificities.append(top_prediction_name)
            except Exception as e:
                print(f"Could not predict A-domain specificities: {e}")
                predicted_specificities = []

        # count number of groups with AMP-domains
        num_amp_binding = len([group for group in groups if "AMP-binding" in group])
        if num_amp_binding != len(predicted_specificities): predicted_specificities = ["amino acid"] * num_amp_binding
        predicted_specificities = iter(predicted_specificities)
        
        # parse the putative domains
        primary_sequence = ["start"]
        for group in groups:
            if "PKS_KS" in group and "AMP-binding" in group: raise ValueError("Found both PKS_KS and AMP-binding in the same domain.")
            elif "PKS_KS" in group:
                polyketide_type = "A"
                if "PKS_ER" in group: polyketide_type = "D"
                elif "PKS_DH" in group: polyketide_type = "C" 
                elif "PKS_KR" in group: polyketide_type = "B"  
                primary_sequence.append(polyketide_type)
            elif "AMP-binding" in group:
                specificity = next(predicted_specificities)
                primary_sequence.append(specificity)
            else:  raise ValueError("Could not find PKS_KS or AMP-binding in the domain.")
        primary_sequence.append("end")

        # return record with parsed primary sequence
        return {   
            "record_id": self.record_id,
            "category": self.category,
            "product": self.product,
            "primary_sequence": primary_sequence,
            "meta_data": {
                "start_proto_cluster": self.start,
                "end_proto_cluster": self.end,
                "product": self.product,
                "category": self.category,
            }
        }
    

class ProtoClusterDict:
    def __init__(self) -> None:
        self._clusters = {}

    def __len__(self) -> int:  
        return len(self._clusters)

    def add_cluster(self, cluster: ProtoCluster) -> None:
        self._clusters[cluster] = cluster

    def assign_gene_to_region(self, gene: Gene) -> None:
        for cluster_hash in self._clusters.values():
            cluster: ProtoCluster = self._clusters[cluster_hash]
            if cluster.start <= gene.start and gene.end <= cluster.end:
                cluster.genes.append(gene)

    @property
    def clusters(self) -> ty.Dict[int, ProtoCluster]:
        return self._clusters
    
    def remove_empty_clusters(self) -> None:
        clusters = self._clusters.copy()
        for cluster in clusters.values():
            if len(cluster) == 0: del self._clusters[cluster]


def get_genes_from_record(record: ty.Dict[str, ty.Any]) -> ty.Dict[str, Gene]:
    """Get genes from an AntiSMASH record."""
    genes = dict()
    features = record["features"]
    for feature in features:
        if feature["type"] != "gene": continue
        location = feature["location"]
        gene_name = feature["qualifiers"].get("gene", None)
        locus_tag = feature["qualifiers"].get("locus_tag", None)
        if locus_tag is not None: name = locus_tag[0]  # prefer locus tag over gene name
        else:
            if gene_name is not None: name = gene_name[0]
            else: raise ValueError("Could not find gene name!")
        if match := re.match(r"\[(\d+):(\d+)\]\(((\+|-))\)", location):
            start = int(match.group(1))
            end = int(match.group(2))
            strand = match.group(3)
        else: continue
        gene = Gene(name=name, start=start, end=end, strand=strand)
        genes[name] = gene
    return genes


def get_protoclusters_from_record(record: ty.Dict[str, ty.Any], categories: ty.List[str] = ["NRPS", "PKS"]) -> ProtoClusterDict:
    protoclusters = ProtoClusterDict()
    for area in record["areas"]:
        for _, protocluster in area["protoclusters"].items():
            category = protocluster["category"]
            product = protocluster["product"]
            if category not in categories: continue
            cluster = ProtoCluster(record_id=record["id"], start=protocluster["core_start"], end=protocluster["core_end"], product=product, category=category)
            protoclusters.add_cluster(cluster)
    return protoclusters


def get_modules_from_record(record: ty.Dict[str, ty.Any]) -> ty.Dict[str, ty.List[str]]:
    """Get modules from an AntiSMASH record."""
    modules_per_gene = defaultdict(list)  # modules ordered in the order they appear in the gene
    modules_type = "antismash.detection.nrps_pks_domains"
    if modules := record["modules"].get(modules_type, None):
        cds_results = modules["cds_results"]
        for gene_name, props in cds_results.items():
            predicted_modules = props["modules"]
            if len(predicted_modules) == 0: continue
            for predicted_module in predicted_modules:
                module_components = predicted_module["components"]
                for module_component in module_components:
                    hit_id = module_component["domain"]["hit_id"]
                    gene_name = module_component["locus"]
                    modules_per_gene[gene_name].append(hit_id)
    return modules_per_gene


def parse_antismash_json(data: ty.Dict[str, ty.Any], predict_specificities: bool = True) -> ty.List[ty.Dict[str, ty.Any]]:
    """Parse AntiSMASH JSON data."""
    records = data["records"]
    queries = []
    for record in records:
        dna = record["seq"]["data"]
        genes = get_genes_from_record(record)
        modules = get_modules_from_record(record)
        protoclusters = get_protoclusters_from_record(record)

        # assign modules to genes
        for gene_name, gene in genes.items(): gene.modules = modules.get(gene_name, [])

        # remove all genes that have no modules
        genes = {name: gene for name, gene in genes.items() if len(gene.modules) > 0}

        # continue if no genes or protoclusters
        if len(protoclusters) == 0 or len(genes) == 0: continue

        # assign genes to clusters
        for gene in genes.values(): protoclusters.assign_gene_to_region(gene)

        # remove clusters without assigned genes
        protoclusters.remove_empty_clusters()

        for protocluster in protoclusters.clusters.values():
            motif_code_query = protocluster.to_motif_code(dna, predict_specificities=predict_specificities)
            queries.append(motif_code_query)

    return queries


def add_mibig(cur, path_mibig_database, path_mibig_jsons):
    """Add MIBiG database to the database."""
    raise NotImplementedError
