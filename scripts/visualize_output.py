import argparse
from retromol.retrosynthesis.result import Result 

parser = argparse.ArgumentParser()
parser.add_argument("--json", required=True, type=str)
args = parser.parse_args()

path = args.json
result = Result.from_json(path)
svg_str = result.draw_molecule(color_monomers=True)

path = path.replace(".json", ".svg")
with open(path, "w", encoding="utf-8") as fo:
    fo.write(svg_str)
