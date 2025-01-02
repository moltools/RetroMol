#!/usr/bin/env python3
import argparse 
import os


def cli():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", "-i", type=str, required=True, help="path to dir with deorphanization results")
    parser.add_argument("--output", "-o", type=str, required=True, help="path to output file")
    return parser.parse_args()


def main():
    args = cli()
    file_paths = [os.path.join(args.input, f) for f in os.listdir(args.input) if f.endswith(".csv")]

    out_file = open(args.output, "w")
    out_file.write("index,file_path\n")

    high_scorers = 0
    records = []
    
    for file_path in file_paths:
        with open(file_path, "r") as f:
            f.readline() # skip header
            _, score, *_ = f.readline().strip().split(",") # only need the first line
            score = float(score)

        # only print 
        if score >= 0.9:
            high_scorers += 1
            # out_file.write(f"{high_scorers},{file_path}\n")
            records.append(file_path)

    # sort file paths and write out to file
    records.sort()
    for i, record in enumerate(records):
        out_file.write(f"{i+1},{record}\n")

    print(f"Number of high scorers: {high_scorers}")
    out_file.close()

if __name__ == "__main__":
    main()
