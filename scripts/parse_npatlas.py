#!/usr/bin/env python3
import argparse
import json
import logging
from collections import Counter

from tqdm import tqdm

#TODO; superclass results and do macrolides

classes_polyketides = [
    "Open-chain polyketides", # INCLUDE
    # "Simple cyclic polyketides", # EXCLUDE
    "Macrolide lactones", # INCLUDE
    "Polyene macrolides", # INCLUDE
    # "Macrolide lactams", # EXCLUDE
    "Erythromycins", # INCLUDE
    "Epothilones", # INCLUDE
    # "Zearalenones", # EXCLUDE
    # "Ansa macrolides",  # EXCLUDE
    "Tylosins", # INCLUDE
    "Bafilomycins", # INCLUDE
    "Antimycins", # INCLUDE
    "Lactam bearing macrolide lactones", # INCLUDE
]

classes_non_ribosomal_peptides = [
    # "Cyclic peptides",
    # "Linear peptides",
    # "Lipopeptides",
    # "Microcystins",
    # "Cyanopeptolins",
    # "Dipeptides",
    # "Tripeptides",
]

def cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str, required=True, help="path to npatlas json file with annotations")
    parser.add_argument("-o", "--output", type=str, required=True, help="path to output csv file")
    parser.add_argument("-l", "--loglevel", type=str, default="INFO", help="log level", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"])
    return parser.parse_args()


def main() -> None:
    args = cli()

    # config logger
    logging.basicConfig(format="[%(asctime)s %(name)s %(levelname)s] %(message)s")
    logger = logging.getLogger(__name__)
    logger.setLevel(args.loglevel)

    # read npatlas json file
    logger.info(f"reading npatlas json file at {args.input}")
    with open(args.input, "r") as f:
        data = json.load(f)
        logger.info(f"loaded {len(data)} records")

    # count instances
    class_counter = Counter()

    # loop over records and extract information for writing to csv
    logger.info(f"writing out to {args.output}")
    out_fo = open(args.output, "w")

    # write header
    out_fo.write("npaid,smiles,classes\n")

    parsed = 0
    for record in tqdm(data, leave=False):
        npaid = record["npaid"]
        smiles = record["smiles"]

        if classifications := record.get("npclassifier", None):
            if class_results := classifications.get("class_results", None):

                # compile desired classes
                desired_classes = classes_polyketides + classes_non_ribosomal_peptides

                # check if class_results overlap with desired_classes, if so write out
                if any(class_result in desired_classes for class_result in class_results):
                    class_results = [class_result for class_result in class_results if class_result in desired_classes]
                    out_fo.write(f"{npaid},{smiles},\"{'|'.join(class_results)}\"\n")
                    parsed += 1

                if class_results is not None:
                    for class_result in class_results:
                        class_counter[class_result] += 1

    out_fo.close()
    logger.info(f"parsed {parsed} records")

    logger.info("done")

    # sort class results and print to stdout
    print("class,count")
    for class_result, count in class_counter.most_common():
        print(f"\"{class_result}\",{count}")

if __name__ == "__main__":
    main()
    