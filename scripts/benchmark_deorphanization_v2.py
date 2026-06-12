import argparse
import glob
import os
import json
from pathlib import Path

from tqdm import tqdm

from retromol_antismash.io import AntiSmashOptions, load_regions
from retromol_antismash.inference.registry import register_domain_model, annotate_region
from retromol_antismash.inference.model_paras import ParasModel
from retromol_antismash.modules import linear_readout
from retromol_database.duckdb import FINGERPRINT_SIZE, open_database
from retromol.model.rules import RuleSet
from retromol_fingerprint.fingerprint import TOKEN_UNK, Vocabulary, Fingerprinter
from retromol_alignment.scoring import create_tanimoto_scoring_matrix, HARDCODED_PK_SCORING
from retromol_alignment.aligner import setup_aligner
from retromol_alignment.pairwise import Converter
from retromol_alignment.ranking import rerank
from retromol_antismash.modules import linear_readout, PKSModule, NRPSModule
from retromol.chem.mol import smiles_to_mol
from retromol.chem.fingerprint import mol_to_morgan_fingerprint, calculate_tanimoto_similarity

import matplotlib.pyplot as plt


# paras predicted amino acids have names indicating stereochemistry, stereochemistry has been removed from these substrates
# in RetroMol's default ruleset
name_dict = {
    "(2S,3R)-2-amino-3-hydroxy-4-(4-nitrophenyl)butanoic acid" : "(2,3)-2-amino-3-hydroxy-4-(4-nitrophenyl)butanoic acid",
    "(2S,6R)-diamino-(5R,7)-dihydroxy-heptanoic acid" : "(2,6)-diamino-(5,7)-dihydroxy-heptanoic acid",
    "(4S)-5,5,5-trichloroleucine" : "(4)-5,5,5-trichloroleucine",
    "(E)-4-methylhex-2-enoic acid" : "4-methylhex-2-enoic acid",
    "(S,E)-2-amino-4-decenoic acid" : "2-amino-4-decenoic acid",
    "2-(1-methylcyclopropyl)-D-glycine" : "2-(1-methylcyclopropyl)-glycine",
    "2R-hydroxy-3-methylpentanoic acid" : "2-hydroxy-3-methylpentanoic acid",
    "2R-hydroxyisovaleric acid" : "2-hydroxyisovaleric acid",
    "2S,3S-diaminobutyric acid" : "2,3-diaminobutyric acid",
    "2S-amino decanoic acid" : "2-amino decanoic acid",
    "2S-amino-4-hexenoic acid" : "2-amino-4-hexenoic acid",
    "2S-amino-8-oxodecanoic acid" : "2-amino-8-oxodecanoic acid",
    "2S-amino-9,10-epoxy-8-oxodecanoic acid" : "2-amino-9,10-epoxy-8-oxodecanoic acid",
    "2S-amino-dodecanoic acid" : "2-amino-dodecanoic acid",
    "2S-amino-octanoic-acid" : "2-amino-octanoic-acid",
    "2S-aminodecanoic acid" : "2-aminodecanoic acid",
    "2S-aminododecanoic acid" : "2-aminododecanoic acid",
    "2S-aminooctanoic acid" : "2-aminooctanoic acid",
    "2S-hydroxyisocaproic acid" : "2-hydroxyisocaproic acid",
    "2S-hydroxyisovaleric acid" : "2-hydroxyisovaleric acid",
    "2S-methyl-3-oxobutyrine" : "2-methyl-3-oxobutyrine",
    "3-methyl-D-aspartic acid wonky" : "3-methylaspartic acid",
    "3R-aminoisobutyric acid" : "3-aminoisobutyric acid",
    "3R-chloroproline" : "3-chloroproline",
    "3R-hydroxy-2,4-diaminobutyric acid" : "3-hydroxy-2,4-diaminobutyric acid",
    "3R-hydroxyasparagine" : "3-hydroxyasparagine",
    "3R-hydroxyaspartic acid" : "3-hydroxyaspartic acid",
    "3R-hydroxyglutamine" : "3-hydroxyglutamine",
    "3R-hydroxyhomotyrosine" : "3-hydroxyhomotyrosine",
    "3R-hydroxyleucine" : "3-hydroxyleucine",
    "3R-methyl-D-aspartic acid wonky" : "3-methylaspartic acid",
    "3R-methylbeta-alanine" : "3-methylbeta-alanine",
    "3R-methylglutamic acid" : "3-methylglutamic acid",
    "3S,4R-dichloroproline" : "3,4-dichloroproline",
    "3S,4S-dihydroxyhomotyrosine" : "3,4-dihydroxyhomotyrosine",
    "3S-aminobutyric acid" : "3-aminobutyric acid",
    "3S-carboxypiperazine" : "3-carboxypiperazine",
    "3S-cyclohex-2-enylalanine" : "3-cyclohex-2-enylalanine",
    "3S-hydroxy-4R-methyloctanoic acid" : "3-hydroxy-4Rmethyloctanoic acid",
    "3S-hydroxy-4S-methylproline" : "3-hydroxy-4-methylproline",
    "3S-hydroxy-6-chlorohistidine" : "3-hydroxy-6-chlorohistidine",
    "3S-hydroxyasparagine" : "3-hydroxyasparagine",
    "3S-hydroxyleucine" : "3-hydroxyleucine",
    "3S-hydroxypipecolic acid" : "3-hydroxypipecolic acid",
    "3S-hydroxyproline" : "3-hydroxyproline",
    "3S-methyl-D-aspartic acid branched" : "3-methylaspartic acid",
    "3S-methyl-D-aspartic acid wonky" : "3-methyl-aspartic acid",
    "3S-methylaspartic acid" : "3-methylaspartic acid",
    "3S-methylaspartic acid branched" : "3-methylaspartic acid",
    "3S-methylleucine" : "3-methylleucine",
    "3S-methylproline" : "3-methylproline",
    "4-hydroxy-D-kynurenine" : "4-hydroxy-kynurenine",
    "4R-E-butenyl-4R-methylthreonine" : "4-butenyl-4-methylthreonine",
    "4R-hydroxyproline" : "4-hydroxyproline",
    "4R-methylproline" : "4-methylproline",
    "4R-propylproline" : "4-propylproline",
    "4S,5-dihydroxy-2S-aminopentanoic acid" : "4,5-dihydroxy-2-aminopentanoic acid",
    "4S-acetyl-5S-methylproline" : "4-acetyl-5-methylproline",
    "4S-hydroxylysine" : "4-hydroxylysine",
    "4S-methylazetidine-2S-carboxylic acid" : "4-methylazetidine-2-carboxylic acid",
    "4S-methylproline" : "4-methylproline",
    "4S-propenylproline" : "4-propenylproline",
    "5S-methylproline" : "5-methylproline",
    "6S-methyl-pipecolic acid" : "6-methyl-pipecolic acid",
    "D-alanine" : "alanine",
    "D-aspartic acid branched" : "aspartic acid",
    "D-glutamic acid branched" : "glutamic acid",
    "D-isovaline" : "isovaline",
    "D-leucine" : "leucine",
    "D-lysergic acid" : "lysergic acid",
    "D-phenylalanine" : "phenylalanine",
    "D-phenyllactic acid" : "phenyllactic acid",
    "D-pipecolic acid" : "pipecolic acid",
    "D-serine" : "serine",
    "R-3-hydroxy-3-methylproline" : "3-hydroxy-3-methylproline",
    "R-aza-beta-tyrosine" : "beta-tyrosine",
    "R-beta-hydroxyphenylalanine" : "beta-hydroxyphenylalanine",
    "R-beta-hydroxytyrosine" : "beta-hydroxytyrosine",
    "R-beta-methylphenylalanine" : "beta-methylphenylalanine",
    "R-beta-methyltryptophan" : "beta-methyltryptophan",
    "R-beta-phenylalanine" : "beta-phenylalanine",
    "R-beta-tyrosine" : "beta-tyrosine",
    "S-adenosylmethionine" : "adenosylmethionine",
    "S-beta-hydroxycyclohex-2S-enylalanine" : "beta-hydroxycyclohex-2-enylalanine",
    "S-beta-hydroxyenduracididine" : "beta-hydroxyenduracididine",
    "S-beta-hydroxyphenylalanine" : "beta-hydroxyphenylalanine",
    "S-beta-methylphenylalanine" : "beta-methylphenylalanine",
    "S-beta-tyrosine" : "beta-tyrosine",

    "3-(3-pyridyl)-alanine": "alanine",
    "2,3-diaminopropionic acid": "2,3-diaminopropionate",
    "O-methyltyrosine" : "N-methyltyrosine",
    "2,4-diaminobutyric acid": "3-hydroxy-2,4-diaminobutyric acid",
    "allo-isoleucine" : "isoleucine",
    "aspartic acid branched" : "aspartic acid",
    "allo-threonine": "threonine",
    "(2,6)-diamino-(5,7)-dihydroxy-heptanoic acid": "2,6-diamino-57-dihydroxy-heptanoic acid",
}


def cli():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mibig-jsons", type=Path, required=True)
    parser.add_argument("--mibig-gbks", type=Path, required=True)
    parser.add_argument("--retromol-db", type=Path, required=True)
    parser.add_argument("--cache-dir", type=Path, required=True)
    return parser.parse_args()



def modules_to_primary_sequence(modules) -> list[list[str]]:
    per_monomer_tokens = []
    for module in modules:
        match module:
            case PKSModule():
                extender_unit_type = module.substrate.extender_unit.value
                per_monomer_tokens.append([extender_unit_type, "PK"])
            case NRPSModule():
                tokens = []
                if substrate := module.substrate:
                    if substrate.name == "unknown":
                        tokens.append(TOKEN_UNK)
                    else:
                        tokens.append(substrate.name)
                else:
                    tokens.append(TOKEN_UNK)
                per_monomer_tokens.append(tokens)
            case _:
                raise ValueError(f"Unknown module {module}")
    return per_monomer_tokens


def main():
    args = cli()

    # Setup fingerprinter
    rules = RuleSet.load_default()
    vocab_tokens = set()
    for rule in tqdm(rules.matching_rules):
        rule_name = rule.name
        vocab_tokens.add(rule_name)
        rule_pseudonyms = rule.pseudonyms
        vocab_tokens.update(set(rule_pseudonyms))
    vocab = Vocabulary(vocab_tokens)
    fingerprinter = Fingerprinter(vocab, n_bits=FINGERPRINT_SIZE, n_hashes=2)

    # Create scoring matrix dynamically, only match on diagonal
    records = []
    for rule in tqdm(rules.matching_rules):
        rule_name = rule.name
        rule_smiles = rule.smiles
        records.append((rule_name, rule_smiles))
    scoring_matrix = create_tanimoto_scoring_matrix(
        records,
        radius=2,
        num_bits=2048,
        stereochemistry=False,
        self_score_tokens=[TOKEN_UNK, "PK_A", "PK_B", "PK_C", "PK_D"],
        self_score=1.0,
        hardcoded_scores=HARDCODED_PK_SCORING,
    )
    print(scoring_matrix.alphabet)
    aligner = setup_aligner(scoring_matrix, mode="global")

    unknowns = []

    def to_identifier(s):
        try:
            return scoring_matrix.alphabet.index(s)
        except:
            unknowns.append(s)
            return scoring_matrix.alphabet.index("<UNK>")

    def from_identifier(i):
        return scoring_matrix.alphabet[i]


    converter = Converter(
        # to_identifier=lambda s: scoring_matrix.alphabet.index(s),
        # from_identifier=lambda i: scoring_matrix.alphabet[i],
        to_identifier=to_identifier,
        from_identifier=from_identifier,
    )

    # Parse mibig accessions with product types + num modules from gbks
    elligible_accessions_gbks = {}
    options = AntiSmashOptions(readout_level="cand_cluster")
    register_domain_model(ParasModel(threshold=0.1, keep_top=1, cache_dir=args.cache_dir))
    for path in tqdm(glob.glob(os.path.join(args.mibig_gbks, "*.gbk"))):
        regions = load_regions(path, options=options)
        if regions:
            region = max(regions, key=lambda r: r.end - r.start)
            annotate_region(region)
            modules = linear_readout(region).modules
            primary_sequence = modules_to_primary_sequence(modules)
            print(primary_sequence)
            num_modules = len(linear_readout(region).modules)
            product_qualifiers = sorted(region.qualifiers.get("product", []))
            if (product_qualifiers == ["NRPS", "T1PKS"] or product_qualifiers == ["NRPS"] or product_qualifiers == ["T1PKS"]) and num_modules >= 2:
                elligible_accessions_gbks[region.id] = primary_sequence
    print(len(elligible_accessions_gbks))

    # Check which accessions have at least one product structure associated with it
    elligible_accessions_jsons = {}
    for path in tqdm(glob.glob(os.path.join(args.mibig_jsons, "*.json"))):
        with open(path, "r") as f:
            data = json.load(f)
        if compounds := data.get("compounds", None):
            for compound in compounds:
                if structure := compound.get("structure", None):
                    elligible_accessions_jsons[data["accession"]] = structure
    print(len(elligible_accessions_jsons))

    # Print len overlap
    overlap = set(elligible_accessions_gbks.keys()) & set(elligible_accessions_jsons.keys())
    print(len(overlap))

    # Setup query objects for deorphanization
    queries = []
    for accession in overlap:
        queries.append((accession, elligible_accessions_gbks[accession], elligible_accessions_jsons[accession]))

    thresholds = [0.4, 0.675]
    ks = [1, 5, 10, 50, 100]
    scores = [[[] for _ in ks] for _ in thresholds]

    with open_database(args.retromol_db, read_only=True) as db:
        for accession, per_monomer_tokens, smiles in tqdm(queries):

            query_mol = smiles_to_mol(smiles)
            query_fp = mol_to_morgan_fingerprint(query_mol, radius=2, num_bits=2048)
            fingerprint = fingerprinter.encode(per_monomer_tokens=per_monomer_tokens)
            primary_sequence = [name_dict.get(x[0], x[0]) for x in per_monomer_tokens]

            print(fingerprint)
            print(primary_sequence)
            print(smiles)

            closest_1000 = db.closest(fingerprint, limit=1000)
            targets = [result.entry.primary_sequence for result in closest_1000]
            reranked = rerank(
                query=primary_sequence,
                targets=targets,
                aligner=aligner,
                converter=converter,
            )
            reranked_results = []
            for result, (score, inverted) in zip(closest_1000, reranked, strict=True):
                reranked_results.append((score, result))
            reranked_results.sort(key=lambda x: x[0], reverse=True)
            reranked_smis = [x[1].entry.raw for x in reranked_results]
            reranked_fps = [mol_to_morgan_fingerprint(smiles_to_mol(smi), radius=2, num_bits=2048) for smi in reranked_smis]
            morgan_scores = [calculate_tanimoto_similarity(query_fp, target_fp) for target_fp in reranked_fps]

            for i, threshold in enumerate(thresholds):
                for j, k in enumerate(ks):
                    top_n = any(score >= threshold for score in morgan_scores[:k])
                    scores[i][j].append(top_n)

    line_plots = [[] for _ in thresholds]
    for i, xs in enumerate(scores):
        for xxs in xs:
            line_plots[i].append(sum(xxs) / len(xxs))

    for line_plot in line_plots:
        print(line_plot)

    print(set(unknowns))

    fig, ax = plt.subplots(figsize=(5, 2))

    ax.plot(
        ks,
        line_plots[0],
        color="#56b4e9",
        marker="o",
        linewidth=2.5,
        markersize=7,
        # label="Tanimoto ≥ 0.4",
        markeredgecolor="black",
        markeredgewidth=1.0,
    )
    ax.plot(
        ks,
        line_plots[1],
        color="#039e73",
        marker="s",
        linewidth=2.5,
        markersize=7,
        # label="Tanimoto ≥ 0.675",
        markeredgecolor="black",
        markeredgewidth=1.0,
    )

    ax.set_xlabel("Top-N")
    ax.set_ylabel("Ratio with a matching compound")
    ax.set_xticks(ks)
    ax.set_xlim(-5, 105)
    ax.set_ylim(0.0, 0.6)

    # ax.legend()
    ax.grid(axis="y", alpha=0.3)

    fig.tight_layout()
    fig.savefig("/Users/davidmeijer/Downloads/top_n_retrieval_ratio.png", dpi=300)
    fig.savefig("/Users/davidmeijer/Downloads/top_n_retrieval_ratio.svg", dpi=300)
    plt.show()


if __name__ == "__main__":
    main()