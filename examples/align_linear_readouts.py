
from retromol.model.submission import Submission
from retromol.model.rules import RuleSet
from retromol.model.readout import LinearReadout
from retromol.pipelines.parsing import run_retromol_with_timeout

# erythromycin
smi1 = r"CC[C@@H]1[C@@]([C@@H]([C@H](C(=O)[C@@H](C[C@@]([C@@H]([C@H]([C@@H]([C@H](C(=O)O1)C)O[C@H]2C[C@@]([C@H]([C@@H](O2)C)O)(C)OC)C)O[C@H]3[C@@H]([C@H](C[C@H](O3)C)N(C)C)O)(C)O)C)C)O)(C)O"
# megalomycin
smi2 = r"CC[C@@H]1[C@@]([C@@H]([C@H](C(=O)[C@@H](C[C@@]([C@@H]([C@H]([C@@H]([C@H](C(=O)O1)C)O[C@H]2C[C@@]([C@H]([C@@H](O2)C)O)(C)O)C)O[C@H]3[C@@H]([C@H](C[C@H](O3)C)N(C)C)O)(C)O[C@H]4C[C@H]([C@H]([C@@H](O4)C)O)N(C)C)C)C)O)(C)O"
# 6-deoxyerythronolide
smi3 = r"CC[C@@H]1[C@@H]([C@@H]([C@H](C(=O)[C@@H](C[C@@H]([C@@H]([C@H]([C@@H]([C@H](C(=O)O1)C)O)C)O)C)C)C)O)C"

smis = [smi1, smi2, smi3]
subs = [Submission(smi) for smi in smis]

ruleset = RuleSet.load_default(match_stereochemistry=False)
results = [run_retromol_with_timeout(sub, ruleset) for sub in subs]
coverages = [res.calculate_coverage() for res in results]
print(coverages)
readouts = [LinearReadout.from_result(res) for res in results]
print(len(readouts))

# TODO:
# align item 2 to item 1: align all extracted paths and pick best, that is first part of alignment
# then take all the unaligned parts and do it again until either has no more paths to align
# parsing results should have few metrics: best scoring part and how long that alignment is, if both are high at least partial good match
# need to sort on this, or have automatic sorting function; need to generalize to other types of linear readouts