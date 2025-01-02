#!/usr/bin/env python3
import argparse 
import json
import os
from tqdm import tqdm


def cli():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str, help="input dir")
    parser.add_argument("-o", "--output", type=str, help="output dir")
    return parser.parse_args()


def main():
    args = cli()

    # get all files in input dir
    input_files = [f for f in os.listdir(args.input) if f.endswith(".json")]

    # loop over input files, check contens, and decide to write out again to out dir
    for input_file in tqdm(input_files):
        with open(os.path.join(args.input, input_file), "r") as file:
            data = json.load(file)

        # only keep records where primary sequence has length 7 or more (includes start and end, so 5 units in total)
        min_target_length = 7
        record_length = len(data["primary_sequence"])
        if record_length < min_target_length:
            continue

        # write out to output dir
        with open(os.path.join(args.output, input_file), "w") as file:
            json.dump(data, file)

    # count files in output dir
    output_files = [f for f in os.listdir(args.output) if f.endswith(".json")]
    print(f"Number of files in output dir: {len(output_files)}")


if __name__ == "__main__":
    main()
