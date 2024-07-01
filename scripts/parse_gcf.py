import argparse
import json
from retromol.npkg.antismash import parse_antismash_json

parser = argparse.ArgumentParser()
parser.add_argument("-i", type=str, required=True, help="path to json")
args = parser.parse_args()

with open(args.i, "r") as f:
    src = json.load(f)

protoclusters = parse_antismash_json(data=src, predict_specificities=False)

print(len(protoclusters))

for i, protocluster in enumerate(protoclusters):
    i += 1
    if i == 16:
        print(protocluster["metaData"]["startProtocluster"])
        print(protocluster["metaData"]["endProtocluster"])
        print(protocluster["recordId"])
