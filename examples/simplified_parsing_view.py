#!/usr/bin/env python3

import argparse

from retromol.chem.mol import encode_mol, smiles_to_mol, smarts_to_mol, mol_to_smiles
from retromol.chem.reaction import smarts_to_reaction
from retromol.model.rules import RuleSet
from retromol.model.submission import Submission
from retromol.model.readout import LinearReadout
from retromol.pipelines.parsing import run_retromol_with_timeout
from retromol.model.readout import MolNode


def cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--smiles", type=str, required=True, help="SMILES string")
    return parser.parse_args()


indicator_smarts = ["[PbH2]", "[SnH2]", "[O]=[C]=[O]"]
indicator_mols = [smarts_to_mol(s) for s in indicator_smarts]
def filter_indicators(products):
    results = []
    for product in products:
        if not any(product.HasSubstructMatch(ind) for ind in indicator_mols):
            results.append(product)
    return results


def main() -> None:
    args = cli()

    # Parse SMILES string with default ruleset
    submission = Submission(args.smiles)
    ruleset = RuleSet.load_default(match_stereochemistry=True)
    result = run_retromol_with_timeout(submission, ruleset)
    coverage = result.calculate_coverage()
    print(f"Compound parsed with coverage {coverage:.2f}")

    rxn_start_capping = smarts_to_reaction("[SnH:1][C:2]~[C:3][PbH:4]>>[OH][C](=O)[C:2]~[C:3][PbH:4].[SnH2:1]")
    rxn_extension_1 = smarts_to_reaction("[C:1][C:2]~[C:3][PbH:4].[SnH:5][C:6]-[C:7][PbH:8]>>[C:1][C:2]~[C:3][C:6]-[C:7][PbH:8].[PbH2:4].[SnH2:5]")
    rxn_extension_2 = smarts_to_reaction("[C:1][C:2]~[C:3][PbH:4].[SnH:5][C:6]=[C:7][PbH:8]>>[C:1][C:2]~[C:3][C:6]=[C:7][PbH:8].[PbH2:4].[SnH2:5]")
    # only works for existing polyketide chain, extender amino acid
    # also need:
    # - existing amino acid chain, extender amino acid
    # - existing amino acid chain, extender polyketide unit
    rxn_extension_3 = smarts_to_reaction("[C:1][C:2]~[C:3][PbH:4].[OH:5][C:6](=[O:7])[C:8][NH2:9]>>[C:1][C:2]~[C:3][C:8][NH2:9].[PbH2:4].[O:5]=[C:6](=[O:7])")
    extension_1_pattern = smarts_to_mol("[SnH][C]-[C][PbH]")
    extension_2_pattern = smarts_to_mol("[SnH][C]=[C][PbH]")
    extension_3_pattern = smarts_to_mol("[OH][C](=O)C[NH2]")
    rxn_end_capping = smarts_to_reaction("[PbH:1][C:2]>>[*][C:2].[PbH2:1]")

    # Retrieve view from result object
    matching_rules = [r for r in ruleset.matching_rules if not r.props.get("polyketide_motif", False)]
    root_enc = encode_mol(submission.mol)
    readout = LinearReadout.from_reaction_graph(root_enc, reaction_graph=result.reaction_graph, identified_only=True)
    building_blocks = []
    for path in readout.paths:
        for item in path:
            name, smiles = item.identity.matched_rule.name, item.identity.matched_rule.smiles if item.identified else (None, None)
            print(name, smiles)
            building_blocks.append(smiles_to_mol(smiles))
    print(building_blocks)

    # cap first building block
    print("START CAPPING")
    products = rxn_start_capping.RunReactants((building_blocks[0],))
    assert len(products) == 1, "Expected exactly 1 product"
    products = list(products[0])
    products = filter_indicators(products)
    assert len(products) == 1, "Expected exactly 1 product"
    product = products[0]

    # TODO: capping start might also be amino acid

    # loop through remaining items and attech to first building block
    for building_block in building_blocks[1:]:
        print("building blocks", mol_to_smiles(product), mol_to_smiles(building_block))
        if building_block.HasSubstructMatch(extension_1_pattern):
            products = rxn_extension_1.RunReactants((product, building_block,))
        elif building_block.HasSubstructMatch(extension_2_pattern):
            products = rxn_extension_2.RunReactants((product, building_block,))
        elif building_block.HasSubstructMatch(extension_3_pattern):
            products = rxn_extension_3.RunReactants((product, building_block,))
        else:
            raise ValueError(f"Unknown building block: {building_block}")
        products = products[0]
        products = filter_indicators(products)
        assert len(products) == 1, "Expected exactly 1 product"
        product = products[0]
        product_smi = mol_to_smiles(product)
        print("result", product_smi, product)

    # end capping
    print("END CAPPING")
    try:
        products = rxn_end_capping.RunReactants((product,))
        assert len(products) == 1, "Expected exactly 1 product"
        products = filter_indicators(products[0])
        product = products[0]
    except:
        product = product
    product_smi = mol_to_smiles(product)
    print("product", product_smi)

    # TODO: can we parse input compound in such a way that it is clear what is start/end for polyketides? want to include start carbonic acid in coverage
    # TODO: allow for hybrids with amino acids
    # TODO: retrieval of sequence of molnodes should break on unknown nodes... and not just paste together identified monomers from assembly graph
    # TODO: use classification for parsing strategy; if predicted to be lipopeptide by NPClassifier then no polyketide rules (add switch to rules)
    # TODO: make interface interactive in app
    # TODO: don't always need capping if amino acid is at start or end..

    # check: gephyronic acid, erythromycin, daptomycin, epothilone
    # check: make sure integration tests for retromol parsing algorithm pass again
    # check: reimplement gui + bionexus into monorepo


if __name__ == "__main__":
    main()
