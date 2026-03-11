import numpy as np
import networkx as nx
import sklearn
from tqdm import tqdm

import matplotlib.pyplot as plt
from sklearn.decomposition import PCA   

from retromol.model.rules import RuleSet
from retromol.model.submission import Submission
from retromol.pipelines.parsing import run_retromol_with_timeout
from retromol_fingerprint.fingerprint import module_graph_fingerprint

compounds = [
    ("erythromycin", r"CC[C@@H]1[C@@]([C@@H]([C@H](C(=O)[C@@H](C[C@@]([C@@H]([C@H]([C@@H]([C@H](C(=O)O1)C)O[C@H]2C[C@@]([C@H]([C@@H](O2)C)O)(C)OC)C)O[C@H]3[C@@H]([C@H](C[C@H](O3)C)N(C)C)O)(C)O)C)C)O)(C)O"),
    ("megalomycin", r"CC[C@@H]1[C@@]([C@@H]([C@H](C(=O)[C@@H](C[C@@]([C@@H]([C@H]([C@@H]([C@H](C(=O)O1)C)O[C@H]2C[C@@]([C@H]([C@@H](O2)C)O)(C)OC(=O)C)C)O[C@@H]3[C@H]([C@@H](C[C@@H](O3)C)N(C)C)O)(C)O[C@H]4C[C@H]([C@H]([C@@H](O4)C)O)N(C)C)C)C)O)(C)O"),
    ("discodermolide", r"C[C@H]1[C@@H](OC(=O)[C@@H]([C@H]1O)C)C[C@@H](/C=C\[C@H](C)[C@@H]([C@@H](C)/C=C(/C)\C[C@H](C)[C@H]([C@H](C)[C@H]([C@@H](C)/C=C\C=C)OC(=O)N)O)O)O"),
    ("dictyostatin", r"C[C@H]1CC[C@H]([C@@H]([C@@H](OC(=O)/C=C\C=C\[C@H]([C@H](C[C@@H](/C=C\[C@@H]([C@@H]([C@H](C1)C)O)C)O)O)C)[C@@H](C)/C=C\C=C)C)O"),
    ("daptomycin", r"CCCCCCCCCC(=O)N[C@@H](CC1=CNC2=CC=CC=C21)C(=O)N[C@H](CC(=O)N)C(=O)N[C@@H](CC(=O)O)C(=O)N[C@H]3[C@H](OC(=O)[C@@H](NC(=O)[C@@H](NC(=O)[C@H](NC(=O)CNC(=O)[C@@H](NC(=O)[C@H](NC(=O)[C@@H](NC(=O)[C@@H](NC(=O)CNC3=O)CCCN)CC(=O)O)C)CC(=O)O)CO)[C@H](C)CC(=O)O)CC(=O)C4=CC=CC=C4N)C"),
    ("daptomycin (no tail)", r"N[C@@H](CC1=CNC2=CC=CC=C21)C(=O)N[C@H](CC(=O)N)C(=O)N[C@@H](CC(=O)O)C(=O)N[C@H]3[C@H](OC(=O)[C@@H](NC(=O)[C@@H](NC(=O)[C@H](NC(=O)CNC(=O)[C@@H](NC(=O)[C@H](NC(=O)[C@@H](NC(=O)[C@@H](NC(=O)CNC3=O)CCCN)CC(=O)O)C)CC(=O)O)CO)[C@H](C)CC(=O)O)CC(=O)C4=CC=CC=C4N)C"),
    ("deoxyerythronolide", r"CC[C@@H]1[C@@H]([C@@H]([C@H](C(=O)[C@@H](C[C@@H]([C@@H]([C@H]([C@@H]([C@H](C(=O)O1)C)O)C)O)C)C)C)O)C")
]

lbls = []
covs = []
grhs = []

ruleset = RuleSet.load_default(match_stereochemistry=False)

found_monomers: set[str] = set()
for name, smi in tqdm(compounds):
    lbls.append(name)

    submission = Submission(smi)
    result = run_retromol_with_timeout(submission, ruleset)

    assembly_graph = result.linear_readout.assembly_graph
    nodes = assembly_graph.monomer_nodes()
    idens = [n.identity.name for n in nodes if n.identified]
    found_monomers.update(idens)

    cov = result.calculate_coverage()
    covs.append(cov)
    
    G = assembly_graph.g   
    for n in G.nodes:
        molnode = G.nodes[n]["molnode"]
        if molnode is not None:
            iden = molnode.identity.name if molnode.identified else "UNIDENTIFIED"
        else:
            iden = "UNIDENTIFIED"
        # Assign iden to node_name_attr
        G.nodes[n]["name"] = iden
    grhs.append(G)

# Create substitution matrix on the fly
print(found_monomers)
found_monomers.add("UNIDENTIFIED")
for x in ("A", "B", "C", "D"):  # these are raw module oxidation states read from a BGC GenBank file
    for i in range(1, 13):
        y = f"{x}{i}"
        found_monomers.add(x)
        found_monomers.add(y)
substitution_matrix = {}
for iden in found_monomers:
    substitution_matrix[iden] = {other_iden: 1.0 if iden == other_iden else 0.0 for other_iden in found_monomers}
substitution_matrix["UNIDENTIFIED"]["UNIDENTIFIED"] = 0.0 
for x in ("A", "B", "C", "D"):
    for i in range(1, 13):
        y = f"{x}{i}"
        substitution_matrix[x][y] = 1.0
        substitution_matrix[y][x] = 1.0
substitution_matrix["B2"]["B6"] = 0.5 
substitution_matrix["B6"]["B2"] = 0.5 
substitution_matrix["D6"]["D2"] = 0.5
substitution_matrix["D2"]["D6"] = 0.5

# Add some synthetic results for parsing GBKs
lbls.append("BGC55")
G = nx.Graph()
G.add_node(0, name="B")
G.add_node(1, name="B")
G.add_node(2, name="A")
G.add_node(3, name="D")
G.add_node(4, name="B")
G.add_node(5, name="B")
G.add_edges_from([(0, 1), (1, 2), (2, 3), (3, 4), (4, 5)])
grhs.append(G)

fprs = []
counted = True
for G in grhs:
    fpr = module_graph_fingerprint(
        graph=G,
        substitution_matrix=substitution_matrix,
        n_bits=2048,
        radius=2,
        node_name_attr="name",
        embedding_dim=16,
        count_similar_edges=True,
        counted=counted,
    )
    fprs.append(fpr)

fprs = np.array(fprs)
print(fprs.shape)

# Center/normalize for PCA
fprs = fprs - np.mean(fprs, axis=0, keepdims=True)
norms = np.linalg.norm(fprs, axis=0, keepdims=True)
fprs = fprs / np.where(norms == 0, 1.0, norms)

# PCA for visualization
pca = PCA(n_components=2)
fprs_2d = pca.fit_transform(fprs)
plt.figure(figsize=(8, 6))
plt.scatter(fprs_2d[:, 0], fprs_2d[:, 1], c='blue')
for i, lbl in enumerate(lbls):
    plt.annotate(lbl, (fprs_2d[i, 0], fprs_2d[i, 1]))
plt.title("PCA of Compound Fingerprints")
plt.xlabel("PC1")
plt.ylabel("PC2")
plt.grid()
plt.show()
