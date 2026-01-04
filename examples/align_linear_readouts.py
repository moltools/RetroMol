
import numpy as np
import matplotlib.pyplot as plt

from retromol.model.submission import Submission
from retromol.model.rules import RuleSet
from retromol.model.readout import LinearReadout
from retromol.pipelines.parsing import run_retromol_with_timeout
from retromol.fingerprint.fingerprint import FingerprintGenerator
from retromol.fingerprint.similarity import calculate_cosine_similarity, calculate_tanimoto_similarity

name1 = "erythromycin"
smi1 = r"CC[C@@H]1[C@@]([C@@H]([C@H](C(=O)[C@@H](C[C@@]([C@@H]([C@H]([C@@H]([C@H](C(=O)O1)C)O[C@H]2C[C@@]([C@H]([C@@H](O2)C)O)(C)OC)C)O[C@H]3[C@@H]([C@H](C[C@H](O3)C)N(C)C)O)(C)O)C)C)O)(C)O"
name2 = "megalomycin"
smi2 = r"CC[C@@H]1[C@@]([C@@H]([C@H](C(=O)[C@@H](C[C@@]([C@@H]([C@H]([C@@H]([C@H](C(=O)O1)C)O[C@H]2C[C@@]([C@H]([C@@H](O2)C)O)(C)O)C)O[C@H]3[C@@H]([C@H](C[C@H](O3)C)N(C)C)O)(C)O[C@H]4C[C@H]([C@H]([C@@H](O4)C)O)N(C)C)C)C)O)(C)O"
name3 = "6-deoxyerythronolide"
smi3 = r"CC[C@@H]1[C@@H]([C@@H]([C@H](C(=O)[C@@H](C[C@@H]([C@@H]([C@H]([C@@H]([C@H](C(=O)O1)C)O)C)O)C)C)C)O)C"
name4 = "gephyronic acid"
smi4 = r"C[C@@H]1[C@@H](O[C@@](C([C@H]1OC)(C)C)([C@@H](C)C[C@H](C)[C@@H]([C@@]2([C@H](O2)[C@@H](C)C=C(C)C)C)O)O)CC(=O)O"
name5 = "nocardichelin B"
smi5 = r"CCCCCCCCCCC/C=C\C(=O)N(CCCCCNC(=O)CCC(=O)N(CCCCCNC(=O)[C@@H]1COC(=N1)C2=CC=CC=C2O)O)O"
name6 = "daptomycin"
smi6 = r"CCCCCCCCCC(=O)N[C@@H](CC1=CNC2=CC=CC=C21)C(=O)N[C@H](CC(=O)N)C(=O)N[C@@H](CC(=O)O)C(=O)N[C@H]3[C@H](OC(=O)[C@@H](NC(=O)[C@@H](NC(=O)[C@H](NC(=O)CNC(=O)[C@@H](NC(=O)[C@H](NC(=O)[C@@H](NC(=O)[C@@H](NC(=O)CNC3=O)CCCN)CC(=O)O)C)CC(=O)O)CO)[C@H](C)CC(=O)O)CC(=O)C4=CC=CC=C4N)C"

names = [name1, name2, name3, name4, name5, name6]
smis = [smi1, smi2, smi3, smi4, smi5, smi6]
subs = [Submission(smi) for smi in smis]

ruleset = RuleSet.load_default(match_stereochemistry=False)
results = [run_retromol_with_timeout(sub, ruleset) for sub in subs]
coverages = [res.calculate_coverage() for res in results]
print(coverages)
readouts = [LinearReadout.from_result(res) for res in results]
print(len(readouts))

generator = FingerprintGenerator(ruleset.matching_rules, tanimoto_threshold=0.6, morgan_radius=2, morgan_num_bits=2048)
print(len(generator.groups), len(generator.monomers))
fingerprints = np.array([generator.fingerprint_from_result(r, num_bits=512, counted=True) for r in results])
print(fingerprints.shape)

def similarity_matrix(fps, sim):
    n = fps.shape[0]
    M = np.zeros((n, n), dtype=float)

    for i in range(n):
        for j in range(n):
            M[i, j] = sim(fps[i], fps[j])
        
    return M

S = similarity_matrix(fingerprints, calculate_cosine_similarity)

plt.imshow(S, vmin=0.0, vmax=1.0, cmap="viridis")
plt.colorbar(label="similarity")
plt.xticks(range(len(names)), names, rotation=45, ha="right")
plt.yticks(range(len(names)), names)
plt.tight_layout()
plt.show()

# TODO:
# align item 2 to item 1: align all extracted paths and pick best, that is first part of alignment
# then take all the unaligned parts and do it again until either has no more paths to align
# parsing results should have few metrics: best scoring part and how long that alignment is, if both are high at least partial good match
# need to sort on this, or have automatic sorting function; need to generalize to other types of linear readouts