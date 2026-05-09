import argparse
import yaml
from pathlib import Path


def cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matching-rules", type=str, help="Path to matching rules file.")
    parser.add_argument("--output", type=str, help="Path to output TSV file.")
    return parser.parse_args()


def main() -> None:
    args = cli()

    output = Path(args.output)

    with open(args.matching_rules, "r") as f_o:
        matching_rules = yaml.safe_load(f_o)

    # Get names for monomers
    names = list(set([r["name"] for r in matching_rules]))
    names = sorted(names)

    with open(output, "w") as f_o:
        header = "\t" + "\t".join(names) + "\n"
        f_o.write(header)

        for name_idx, name in enumerate(names):
            scores = [0.0] * len(names)

            # Fill in self-score
            scores[name_idx] = 1.0

            # Write out scores
            f_o.write(name + "\t" + "\t".join(map(str, scores)) + "\n")


if __name__ == "__main__":
    main()
