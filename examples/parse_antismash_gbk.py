#!/usr/bin/env python3

import argparse
import glob
import os
from pathlib import Path

from retromol_antismash.io import AntiSmashOptions, load_regions
from retromol_antismash.inference.registry import (
    DOMAIN_MODELS,
    GENE_MODELS,
    register_domain_model,
    register_gene_model,
    annotate_region,
)
from retromol_antismash.inference.model_paras import ParasModel
from retromol_antismash.inference.model_pfam import PfamModel
from retromol_antismash.modules import NRPSModule, PKSModule, linear_readout


def cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gbk", required=True, type=str, help="Path to antiSMASH cluster GBK file.")
    parser.add_argument("--hmms", nargs="+", required=False, type=str, default=None, help="Path to HMM Pfam cluster files.")
    return parser.parse_args()


def pprint_module(module: PKSModule | NRPSModule) -> str:
    """
    Pretty print a module.

    :param module: Module to print to stdout.
    :return: String representation of the module.
    """
    if isinstance(module, PKSModule):
        return f"PKS Module: extender_unit={module.substrate.extender_unit.value}, AT_loading_mode={module.anatomy.AT_loading_mode}"
    elif isinstance(module, NRPSModule):
        return f"NRPS Module: substrate={module.substrate.name}"
    else:
        return "Unknown Module"


def main() -> None:
    args = cli()

    # Register domain models
    register_domain_model(ParasModel(threshold=0.1, keep_top=3))

    # Register gene models
    if args.hmms:
        hmm_files = glob.glob(os.path.join(args.hmms, "*.hmm"))
        for hmm_file in hmm_files:
            label = Path(os.path.basename(hmm_file)).stem
            register_gene_model(PfamModel(hmm_path=hmm_file, label=label))

    print(f"registered domain models: {list(DOMAIN_MODELS)}")
    print(f"registered gene models: {list(GENE_MODELS)}")

    # Parse regions from GBK file
    options = AntiSmashOptions()
    regions = load_regions(args.gbk, options)
    print(f"Found {len(regions)} regions in {args.gbk}")

    # Annotate regions using the domain and gene models
    by_orf = True
    for region in regions:
        annotate_region(region)

        # Output some information about the parsed regions
        for r in regions:
            print("region id:", region.id, sep="\t")
            print("file name:", region.file_name, sep="\t")
            readout = linear_readout(r)
            print(readout)
            print(readout.modifiers)
            result = readout.biosynthetic_order(by_orf=by_orf)
            if by_orf:
                for gene_id, modules in result:
                    for mi, module in enumerate(modules):
                        print(f"{gene_id}\t{mi + 1}\t{pprint_module(module)}")
            else:
                for mi, module in enumerate(result):
                    print(f"{mi + 1}\t{pprint_module(module)}")


if __name__ == "__main__":
    main()
