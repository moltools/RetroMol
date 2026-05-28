import argparse

import numpy as np

from retromol_antismash.io import AntiSmashOptions, load_regions
from retromol_antismash.inference.registry import register_domain_model, annotate_region
from retromol_antismash.inference.model_paras import ParasModel
from retromol_antismash.modules import linear_readout


def cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gbk", required=True)
    parser.add_argument("--cache-dir", required=True)
    return parser.parse_args()


def main() -> None:
    args = cli()

    options = AntiSmashOptions(readout_level="cand_cluster")
    register_domain_model(ParasModel(threshold=0.1, keep_top=1, cache_dir=args.cache_dir))
    regions = load_regions(args.gbk, options=options)
    lengths = [r.end - r.start for r in regions]
    longest_region_idx = np.argmax(lengths)
    region = regions[longest_region_idx]
    annotate_region(region)
    modules = linear_readout(region).modules
    modules = sorted(modules, key=lambda x: x.start, reverse=False)

    for module in modules:
        print(module.start, module.end, module.substrate)


if __name__ == "__main__":
    main()
