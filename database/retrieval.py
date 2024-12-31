#!/usr/bin/env python3
import argparse
import typing as ty
import os
import re
import json
import joblib

import pandas as pd
from scipy.stats import chi2_contingency
from statsmodels.stats.multitest import multipletests
from collections import defaultdict
from Bio.Seq import Seq  # biopython
from parasect.api import run_parasect # pip install git+https://github.com/BTheDragonMaster/parasect.git@webapp
from tqdm import tqdm
from versalign.motif import Motif
from versalign.pairwise import PairwiseAlignment, align_pairwise
from versalign.sequence import Sequence
from populate_database import connect_to_database


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
                        
                        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                        # TODO: make sure parsing is better:
                        # improvement 1: get AA sequence directly from features list
                        # improvement 2: sometimes number of predicted domains does not align with number of units we are expecting
                        #                to mitigate this group number of expected domains per gene
                        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                         
                        gene_start = gene.start
                        gene_end = gene.end
                        gene_strand = gene.strand
                        gene_dna = dna[gene_start:gene_end]
                        if gene_strand == "-": gene_dna = str(Seq(gene_dna).reverse_complement())
                        # check if seq is divisible by 3, otherwise add trailing Ns
                        if len(gene_dna) % 3 != 0: gene_dna += "N" * (3 - len(gene_dna) % 3)
                        gene_protein = str(Seq(gene_dna).translate())
                        gene_fasta = f">{gene.name}\n{gene_protein}"
                        to_predict.append(gene_fasta)
                        gene_sequence_added = True

                        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
                for j, result in enumerate(results):
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
            "a_domain_fasta": to_predict,
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

        # get taxonomy data
        annotations = record["annotations"]
        source = annotations["source"]
        organism = annotations["organism"]
        taxonomy = annotations["taxonomy"]

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
            motif_code_query["source"] = source
            motif_code_query["organism"] = organism
            motif_code_query["taxonomy"] = taxonomy
            queries.append(motif_code_query)

    return queries


def parse_scoring_matrix(path_to_file):
    scoring_matrix = {}
    # top row and first column are both headers, names might contain commas as well
    # open with pandas and parse from there to avoid issues with commas in names
    df = pd.read_csv(path_to_file, index_col=0)
    for name in df.columns:
        for other_name, score in df[name].items():
            scoring_matrix[(name, other_name)] = score
    return scoring_matrix


def cli() -> argparse.Namespace:
    """Command line interface."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--scoring-matrix", type=str, required=True, help="path to scoring matrix csv file")
    parser.add_argument("--antismash", type=str, required=True, help="path to antismash output json file")
    return parser.parse_args()


class SequenceMotif(Motif):
    def __init__(self, name: str) -> None: self._name = name
    def __eq__(self, other):
        if not isinstance(other, SequenceMotif): return False
        return self._name == other._name
    def __str__(self): return f"{self._name}"
    def __hash__(self): return hash(self._name)
    @property
    def name(self) -> str:return self._name


def get_targets(organism_types, organism_genus) -> ty.List[Sequence]:
    """Get targets from the database."""
    cur, conn = connect_to_database()
    query = """
        SELECT STRING_AGG(m.motif_name, '|' ORDER BY psm.position) AS primary_sequence, STRING_AGG(DISTINCT cps.compound_id, ',') AS compound_ids 
        FROM primary_sequences_motifs psm JOIN motifs m ON psm.motif_id = m.motif_name 
        JOIN compounds_primary_sequences cps ON psm.primary_sequence_id = cps.primary_sequence_id 
        JOIN compounds_organisms co ON cps.compound_id = co.compound_id JOIN organisms o ON co.organism_id = o.id 
        LEFT JOIN protoclusters_primary_sequences pps ON psm.primary_sequence_id = pps.primary_sequence_id 
        WHERE pps.primary_sequence_id IS NULL AND (%s IS NULL OR o.type = ANY(%s)) AND (%s IS NULL OR o.genus = ANY(%s)) 
        GROUP BY psm.primary_sequence_id;
    """
    cur.execute(query, (organism_types, organism_types, organism_genus, organism_genus))
    results = cur.fetchall()
    cur.close()
    conn.close()
    sequences = []
    primary_sequence_to_compound_ids = {}
    for i, (primary_sequence, compound_ids) in enumerate(results):
        sequence = Sequence(identifier=f"target_{i}", motifs=[SequenceMotif(name) for name in primary_sequence.split("|")])
        sequences.append(sequence)
        primary_sequence_to_compound_ids[primary_sequence] = compound_ids.split(",")
    return sequences, primary_sequence_to_compound_ids


def enrich_for_bioactivity(compound_ids: ty.Set[str]) -> None:
    """Enrich the selected compounds for bioactivity."""
    compound_ids = list(compound_ids)
    cur, conn = connect_to_database()
    # Query to get bioactivity counts for compounds in list and not in list
    counts_query = """
        SELECT ba.bioactivity_type, 
               COUNT(*) FILTER (WHERE ba.compound_id = ANY(%s)) AS count_in_list, 
               COUNT(*) FILTER (WHERE ba.compound_id != ALL(%s)) AS count_not_in_list 
        FROM compounds_bioactivities ba 
        GROUP BY ba.bioactivity_type;
    """
    # Query to get total counts of bioactivity associations in list and not in list
    total_query = """
        SELECT 
            COUNT(*) FILTER (WHERE compound_id = ANY(%s)) AS total_in_list,
            COUNT(*) FILTER (WHERE compound_id != ALL(%s)) AS total_not_in_list
        FROM compounds_bioactivities;
    """
    # Execute counts query
    cur.execute(counts_query, (compound_ids, compound_ids))
    counts = cur.fetchall()
    
    # Execute total counts query
    cur.execute(total_query, (compound_ids, compound_ids))
    totals = cur.fetchone()
    total_in_list = totals[0]
    total_not_in_list = totals[1]

    cur.close()
    conn.close()

    # Perform chi-squared test for each bioactivity_type
    enrichment_results = []
    p_values = []
    bioactivity_types = []
    for bioactivity_type, count_in, count_not_in in counts:
        a = count_in
        b = len(compound_ids) - a
        c = total_in_list - a
        d = total_not_in_list - b
        contingency_table = [[a, b], [c, d]]
        try:
            if count_in > 0: 
                chi2, p, dof, ex = chi2_contingency(contingency_table)
                enrichment_results.append({
                    'bioactivity_type': bioactivity_type,
                    'count_in_list': a,
                    'count_not_in_list': b,
                    'p_value': p
                })
                p_values.append(p)
                bioactivity_types.append(bioactivity_type)
        except Exception as e:
            print(f"Could not perform chi-squared test for bioactivity type {bioactivity_type}: {e}")
    
    if len(p_values) > 0:
        # Apply Bonferroni-Hochberg correction
        reject, pvals_corrected, _, _ = multipletests(p_values, method='fdr_bh')

        # Update results with corrected p-values
        for i in range(len(enrichment_results)):
            enrichment_results[i]['p_value_corrected'] = pvals_corrected[i]
            enrichment_results[i]['reject_null'] = reject[i]

        print(f"Enrichment results for bioactivity:")
        for result in enrichment_results:
            print(f"Bioactivity type: {result['bioactivity_type']}, count in list: {result['count_in_list']}, count not in list: {result['count_not_in_list']}, p-value: {result['p_value']}")


def main() -> None:
    """Main function."""
    args = cli()
    scoring_matrix = parse_scoring_matrix(args.scoring_matrix)
    def score_func(a, b) -> int: return scoring_matrix[(str(a.name), str(b.name))]
    with open(args.antismash, "r") as file: data = json.load(file)
    queries = parse_antismash_json(data, predict_specificities=True)
    print(f"number of queries: {len(queries)}") 
    organism_types = ["bacterium"]
    organism_genus = ["Streptomyces"]
    print(f"retrieving targets for organism types {organism_types} and organism genus {organism_genus}")
    targets, primary_sequence_to_compound_ids = get_targets(organism_types if len(organism_types) else None, organism_genus if len(organism_genus) else None)
    print(f"number of database targets: {len(targets)}")
    for query in queries:
        record_id = query["record_id"]
        primary_sequence = query["primary_sequence"]
        print(f"query primary sequence: {primary_sequence}")
        query = Sequence(identifier=record_id, motifs=[SequenceMotif(name) for name in primary_sequence])
        scores = []
        for target in tqdm(targets):
            _, _, score = align_pairwise(
                seq_a=query, 
                seq_b=target,
                score_func=score_func,
                algorithm=PairwiseAlignment.SMITH_WATERMAN,
                options={"gap_penalty": 3}
            )
            scores.append((target, score))
        scores = sorted(scores, key=lambda x: x[1], reverse=True)
        colleced_compound_ids = set()
        for i, (target, score) in enumerate(scores[:50]):
            target_primary_sequence = "|".join([str(m) for m in target._motifs])
            target_compound_ids = primary_sequence_to_compound_ids.get(target_primary_sequence, [])
            colleced_compound_ids.update(target_compound_ids)
            print(f"({i}) record: {record_id}, target: {target.identifier}, score: {score}, compound_ids: {target_compound_ids}")
        print(f"collected compound ids: {colleced_compound_ids}")
        enrich_for_bioactivity(colleced_compound_ids)
    exit(1)


if __name__ == "__main__":
    main()
