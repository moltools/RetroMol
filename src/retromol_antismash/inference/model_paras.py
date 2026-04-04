"""Module for the PARAS domain inference model."""

import logging
import os
from dataclasses import dataclass
from importlib.resources import files
from pathlib import Path
from typing import Any, TYPE_CHECKING

import joblib
import numpy as np

import retromol_antismash.data
from retromol_antismash.inference.base import DomainInferenceModel
from retromol_antismash.model import Domain, InferenceResult
from retromol_antismash.utils import download_and_prepare

try:
    from pyhmmer import easel, plan7, hmmer
    PYHMMER_AVAILABLE = True

    HMM_DB_PATH = str(files(retromol_antismash.data).joinpath("AMP-binding_converted.hmm"))
    with plan7.HMMFile(HMM_DB_PATH) as hmm_file:
        HMM_DB = list(hmm_file)

except ImportError:
    plan7 = None
    PYHMMER_AVAILABLE = False
    HMM_DB = None

if TYPE_CHECKING:
    from pyhmmer import easel, plan7, hmmer

log = logging.getLogger(__name__)


PARAS_CACHE_DIR = os.getenv("PARAS_CACHE_DIR", "paras_cache")
PARAS_DOWNLOAD_URL = "https://zenodo.org/records/17224548/files/all_substrates_model.paras.gz?download=1"


_PARAS_MODEL_PATH_CACHE: dict[str, Path] = {}
_PARAS_MODEL_OBJ_CACHE: dict[str, object] = {}

VALID = set("ACDEFGHIKLMNPQRSTVWY-")


FEATURE_NAMES = [
    "WOLS870101",
    "WOLS870102",
    "WOLS870103",
    "FAUJ880109",
    "GRAR740102",
    "RADA880108",
    "ZIMJ680103",
    "TSAJ990101",
    "CHOP780201",
    "CHOP780202",
    "CHOP780203",
    "ZIMJ680104",
    "NEU1",
    "NEU2",
    "NEU3",
]


FEATURES = {
    "-": [0.00, 0.00, 0.00, 1, 8.3, 0.21, 13.59, 145.2, 1.00, 1.03, 0.99, 6.03, 0.06, 0.00, 0.10], 
    "A": [0.07, -1.73, 0.09, 0, 8.1, -0.06, 0.00, 90.0, 1.42, 0.83, 0.66, 6.00, 0.06, -0.25, 0.25], 
    "C": [0.71, -0.97, 4.13, 0, 5.5, 1.36, 1.48, 103.3, 0.70, 1.19, 1.19, 5.05, -0.56, -0.40, -0.14], 
    "D": [3.64, 1.13, 2.36, 1, 13.0, -0.80, 49.70, 117.3, 1.01, 0.54, 1.46, 2.77, 0.97, -0.08, 0.08], 
    "E": [3.08, 0.39, -0.07, 1, 12.3, -0.77, 49.90, 142.2, 1.51, 0.37, 0.74, 3.22, 0.85, -0.10, -0.05], 
    "F": [-4.92, 1.30, 0.45, 0, 5.2, 1.27, 0.35, 191.9, 1.13, 1.38, 0.60, 5.48, -0.99, 0.18, 0.15], 
    "G": [2.23, -5.36, 0.30, 0, 9.0, -0.41, 0.00, 64.9, 0.57, 0.75, 1.56, 5.97, 0.32, -0.32, 0.28], 
    "H": [2.41, 1.74, 1.11, 1, 10.4, 0.49, 51.60, 160.0, 1.00, 0.87, 0.95, 7.59, 0.15, -0.03, -0.10], 
    "I": [-4.44, -1.68, -1.03, 0, 5.2, 1.31, 0.13, 163.9, 1.08, 1.60, 0.47, 6.02, -1.00, -0.03, 0.10], 
    "K": [2.84, 1.41, -3.14, 2, 11.3, -1.18, 49.50, 167.3, 1.16, 0.74, 1.01, 9.74, 1.00, 0.32, 0.11], 
    "L": [-4.19, -1.03, -0.98, 0, 4.9, 1.21, 0.13, 164.0, 1.21, 1.30, 0.59, 5.98, -0.83, 0.05, 0.01], 
    "M": [-2.49, -0.27, -0.41, 0, 5.7, 1.27, 1.43, 167.0, 1.45, 1.05, 0.60, 5.74, -0.68, -0.01, 0.04], 
    "N": [3.22, 1.45, 0.84, 2, 11.6, -0.48, 3.38, 124.7, 0.67, 0.89, 1.56, 5.41, 0.70, -0.06, 0.17], 
    "P": [-1.22, 0.88, 2.23, 0, 8.0, 1.1, 1.58, 122.9, 0.57, 0.55, 1.52, 6.30, 0.45, 0.23, 0.41], 
    "Q": [2.18, 0.53, -1.14, 2, 10.5, -0.73, 3.53, 149.4, 1.11, 1.10, 0.98, 5.65, 0.71, -0.02, 0.12], 
    "R": [2.88, 2.52, -3.44, 4, 10.5, -0.84, 52.00, 194.0, 0.98, 0.93, 0.95, 10.76, 0.80, 0.19, -0.41], 
    "S": [1.96, -1.63, 0.57, 1, 9.2, -0.50, 1.67, 95.4, 0.77, 0.75, 1.43, 5.68, 0.48, -0.15, 0.23], 
    "T": [0.92, -2.09, -1.40, 1, 8.6, -0.27, 1.66, 121.5, 0.83, 1.19, 0.96, 5.66, 0.38, -0.10, 0.29], 
    "V": [-2.69, -2.53, -1.29, 0, 5.9, 1.09, 0.13, 139.0, 1.06, 1.70, 0.50, 5.96, -0.75, -0.19, 0.03], 
    "W": [-4.75, 3.65, 0.85, 1, 5.4, 0.88, 2.10, 228.2, 1.08, 1.37, 0.96, 5.89, -0.57, 0.31, 0.34], 
    "Y": [1.39, 2.32, 0.01, 1, 6.2, 0.33, 1.61, 197.0, 0.69, 1.47, 1.14, 5.66, -0.35, 0.40, -0.02], 
}

POSITIONS_ACTIVE_SITE = [
    13,
    16,
    17,
    41,
    46,
    47,
    48,
    49,
    50,
    51,
    52,
    55,
    93,
    94,
    125,
    126,
    127,
    128,
    129,
    152,
    153,
    154,
    155,
    156,
    157,
    158,
    159,
    160,
    161,
    162,
    163,
    164,
    165,
    166,
]


LABEL_TO_SMILES = {
    "(2S,3R)-2-amino-3-hydroxy-4-(4-nitrophenyl)butanoic acid": r"C1=CC(=CC=C1C[C@H]([C@@H](C(=O)O)N)O)[N+](=O)[O-]",
    "(2S,6R)-diamino-(5R,7)-dihydroxy-heptanoic acid": r"C(C[C@@H](C(=O)O)N)[C@H]([C@@H](CO)N)O",
    "(4S)-5,5,5-trichloroleucine": r"CCC(=O)CCCCC[C@@H](C(=O)O)N",
    "(E)-4-methylhex-2-enoic acid": r"CCC(C)/C=C/C(=O)O",
    "(S,E)-2-amino-4-decenoic acid": r"CCCCC/C=C/C[C@@H](C(=O)O)N",
    "1-(1,1-dimethylallyl)-tryptophan": r"CC(C)(C=C)N1C=C(C2=CC=CC=C21)C[C@@H](C(=O)O)N",
    "1-aminocyclopropane-1-carboxylic acid": r"C(O)(=O)C1(CC1)(N)",
    "1-pyrroline-5-carboxylic acid": r"O=C(O)C1/N=C\CC1",
    "10,14-dimethyloctadecanoic acid": r"OC(CCCCCCCCC(C)CCCC(C)CCCC)=O",
    "2,3-diaminobutyric acid": r"NC(C)[C@@H](C(=O)O)N",
    "2,3-diaminopropionic acid": r"C([C@@H](C(=O)O)N)N",
    "2,3-dihydroxy-para-aminobenzoic acid": r"C1=CC(=C(C(=C(N)1)O)O)C(=O)O",
    "2,3-dihydroxybenzoic acid": r"C1=CC(=C(C(=C1)O)O)C(=O)O",
    "2,3-dihydroxyhexadecanoic acid": r"CCCCCCCCCCCCCC(C(C(=O)O)O)O",
    "2,4-diaminobutyric acid": r"C(CN)[C@@H](C(=O)O)N",
    "2,4-dihydroxypentanoic acid": r"CC(CC(C(=O)O)O)O",
    "2-(1-methylcyclopropyl)-D-glycine": r"CC1(CC1)[C@H](C(=O)O)N",
    "2-amino-3,5-dimethyl-4-hexenoic Acid": r"CC(C=C(C)C)C(C(=O)O)N",
    "2-amino-3-hydroxycyclopent-2-enone": r"C1CC(=O)C(=C1O)N",
    "2-amino-6-hydroxy-4-methyl-8-oxodecanoic acid": r"CCC(=O)CC(CC(C)CC(C(=O)O)N)O",
    "2-aminoadipic acid": r"C(C[C@@H](C(=O)O)N)CC(=O)O",
    "2-aminobutyric acid": r"CC[C@@H](C(=O)O)N",
    "2-aminoisobutyric acid": r"O=C(O)C(N)(C)C",
    "2-carboxy-6-hydroxyoctahydroindole": r"N1[C@H](C(=O)O)C[C@@H]2CC[C@@H](O)C[C@H]12",
    "2-chloro-3,5-dihydroxy-4-methylphenylglycine": r"CC1=C(O)C(Cl)=C(C=C(O)1)[C@@H](C(=O)O)N",
    "2-chlorobenzoic acid": r"C1=CC=C(C(=C1)C(=O)O)Cl",
    "2-hydroxy-4-methylpentanoic acid": r"CC(C)CC(C(=O)O)O",
    "2-hydroxypent-4-enoic acid": r"C=CCC(C(=O)O)O",
    "2-ketoglutaric acid": r"C(CC(=O)O)C(=O)C(=O)O",
    "2-ketoisocaproic acid": r"O=C(C(=O)O)CC(C)C",
    "2-ketoisovaleric acid": r"O=C(C(=O)O)C(C)C",
    "2-methylserine": r"C[C@](CO)(C(=O)O)N",
    "2-sulfamoylacetic acid": r"C(C(=O)O)S(=O)(=O)N",
    "2R-hydroxy-3-methylpentanoic acid": r"CCC(C)[C@H](C(=O)O)O",
    "2R-hydroxyisovaleric acid": r"CC(C)[C@H](C(=O)O)O",
    "2S,3S-diaminobutyric acid": r"C[C@@H]([C@@H](C(=O)O)N)N",
    "2S-amino decanoic acid": r"CCCCCCCC[C@H](N)C(=O)O",
    "2S-amino-4-hexenoic acid": r"C/C=C/CC(C(=O)O)N",
    "2S-amino-8-oxodecanoic acid": r"CCC(=O)CCCCC[C@@H](C(=O)O)N",
    "2S-amino-9,10-epoxy-8-oxodecanoic acid": r"C1C(O1)C(=O)CCCCC[C@@H](C(=O)O)N",
    "2S-amino-dodecanoic acid": r"CCCCCCCCCC[C@@H](C(=O)O)N",
    "2S-amino-octanoic-acid": r"CCCCCC[C@@H](C(=O)O)N",
    "2S-aminodecanoic acid": r"CCCCCCCC[C@@H](C(=O)O)N",
    "2S-aminododecanoic acid": r"CCCCCCCCCC[C@@H](C(=O)O)N",
    "2S-aminooctanoic acid": r"CCCCCC[C@@H](C(=O)O)N",
    "2S-hydroxyisocaproic acid": r"CC(C)C[C@@H](C(=O)O)O",
    "2S-hydroxyisovaleric acid": r"CC(C)[C@@H](C(=O)O)O",
    "2S-methyl-3-oxobutyrine": r"CC(=O)[C@](C)(N)C(=O)O",
    "3,3-dihomo-4-methoxytyrosine": r"N[C@@H](CCCC1=CC=C(OC)C=C1)C(=O)O",
    "3,3-dihomophenylalanine": r"N[C@@H](CCCC1=CC=CC=C1)C(=O)O",
    "3,3-dihomotyrosine": r"N[C@@H](CCCC1=CC=C(O)C=C1)C(=O)O",
    "3,4-dehydrolysine": r"C(CCN)=C[C@@H](C(=O)O)N",
    "3,4-dihydroxybenzoic acid": r"C1=CC(=C(C=C1C(=O)O)O)O",
    "3,5-dichloro-4-hydroxyphenylglycine": r"C1=C(Cl)C(=C(Cl)C=C1[C@@H](C(=O)O)N)O",
    "3,5-dihydroxyphenylglycine": r"N[C@H](C(=O)O)c1cc(O)cc(O)c1",
    "3-(2-nitrocyclopropylalanine)": r"C1[C@H]([C@@H]1[N+](=O)[O-])C[C@@H](C(=O)O)N",
    "3-(3-pyridyl)-alanine": r"C1=CC(=CN=C1)C[C@@H](C(=O)O)N",
    "3-amino-2,4-dihydroxybenzoic acid": r"C1=CC(=C(C(=C1C(=O)O)O)N)O",
    "3-amino-4-hydroxybenzoic acid": r"C1=CC(=C(C=C1C(=O)O)N)O",
    "3-amino-6-hydroxy-2-piperidone": r"C1CC(NC(=O)C1N)O",
    "3-aminoisobutyric acid": r"CC(CN)C(=O)O",
    "3-chlorotyrosine": r"C1=C(Cl)C(=CC=C1C[C@@H](C(=O)O)N)O",
    "3-hydroxy-4-methylproline": r"CC1C(O)[C@H](NC1)C(=O)O",
    "3-hydroxy-O-methyl-5-methyltyrosine": r"C1=C(O)C(=C(C)C=C1C[C@@H](C(=O)O)N)OC",
    "3-hydroxy-O-methyltyrosine": r"C1=C(O)C(=CC=C1C[C@@H](C(=O)O)N)OC",
    "3-hydroxy-para-aminobenzoic acid": r"C1=CC(=C(C=C1C(=O)O)O)N",
    "3-hydroxyasparagine": r"N[C@H](C(O)=O)C(O)C(N)=O",
    "3-hydroxyaspartic acid": r"N[C@@H](C(C(=O)O)O)(C(=O)O)",
    "3-hydroxyglutamine": r"C(C([C@@H](C(=O)O)N)O)C(=O)N",
    "3-hydroxykynurenine": r"C1=CC(=C(C(=C1)O)N)C(=O)C[C@@H](C(=O)O)N",
    "3-hydroxyleucine": r"CC(C)C([C@@H](C(=O)O)N)O",
    "3-hydroxypicolinic acid": r"C1=CC(=C(N=C1)C(=O)O)O",
    "3-hydroxyquinaldic acid": r"c1ccc2c(c1)cc(c(n2)C(=O)O)O",
    "3-hydroxytyrosine": r"C1=CC(=C(C=C1C[C@@H](C(=O)O)N)O)O",
    "3-hydroxyvaline": r"CC(O)(C)[C@@H](C(=O)O)N",
    "3-methoxyanthranilic acid": r"COC1=CC=CC(=C1N)C(=O)O",
    "3-methoxyaspartic acid": r"N[C@H](C(C(=O)O)OC)(C(=O)O)",
    "3-methyl-D-aspartic acid wonky": r"C[C@@H]([C@H](C(=O)O)N)C(=O)O",
    "3-methylasparagine": r"CC([C@@H](C(=O)O)N)C(=O)N",
    "3-methylaspartic acid": r"CC([C@@H](C(=O)O)N)C(=O)O",
    "3-methylglutamic acid": r"CC(CC(=O)O)[C@@H](C(=O)O)N",
    "3-methylleucine": r"CC(C)C(C)[C@@H](C(=O)O)N",
    "3-nitrotyrosine": r"C1=CC(=C(C=C1C[C@@H](C(=O)O)N)[N+](=O)[O-])O",
    "3R-aminoisobutyric acid": r"C[C@H](CN)C(=O)O",
    "3R-chloroproline": r"C1[C@@H](Cl)[C@H](NC1)C(=O)O",
    "3R-hydroxy-2,4-diaminobutyric acid": r"NC[C@@H](O)[C@@H](C(=O)O)N",
    "3R-hydroxyasparagine": r"N[C@H](C(O)=O)[C@@H](O)C(N)=O",
    "3R-hydroxyaspartic acid": r"N[C@@H]([C@H](C(=O)O)O)(C(=O)O)",
    "3R-hydroxyglutamine": r"C([C@H]([C@@H](C(=O)O)N)O)C(=O)N",
    "3R-hydroxyhomotyrosine": r"C1=CC(=CC=C1C[C@H]([C@@H](C(=O)O)N)O)O",
    "3R-hydroxyleucine": r"CC(C)[C@H]([C@@H](C(=O)O)N)O",
    "3R-methyl-D-aspartic acid wonky": r"N[C@H]([C@@H](O)C(O)=O)C(O)=O",
    "3R-methylbeta-alanine": r"NC[C@@H](C)C(=O)O",
    "3R-methylglutamic acid": r"C[C@H](CC(=O)O)[C@@H](C(=O)O)N",
    "3S,4R-dichloroproline": r"Cl[C@H]1[C@@H](Cl)[C@H](NC1)C(=O)O",
    "3S,4S-dihydroxyhomotyrosine": r"C1=CC(=CC=C1[C@H](O)[C@H]([C@@H](C(=O)O)N)O)O",
    "3S-aminobutyric acid": r"C[C@@H](CC(=O)O)N",
    "3S-carboxypiperazine": r"C1NN[C@H](C(=O)O)CC1",
    "3S-cyclohex-2-enylalanine": r"C1C=C[C@H](CC1)C[C@@H](C(=O)O)N",
    "3S-hydroxy-4R-methyloctanoic acid": r"CCCC[C@H]([C@H](CC(O)=O)O)C",
    "3S-hydroxy-4S-methylproline": r"C[C@@H]1[C@H](O)[C@H](NC1)C(=O)O",
    "3S-hydroxy-6-chlorohistidine": r"C1=C(NC(Cl)=N1)[C@H]([C@@H](C(=O)O)N)O",
    "3S-hydroxyasparagine": r"N[C@H](C(O)=O)[C@H](O)C(N)=O",
    "3S-hydroxyleucine": r"CC(C)[C@@H]([C@@H](C(=O)O)N)O",
    "3S-hydroxypipecolic acid": r"C1C[C@@H]([C@H](NC1)C(=O)O)O",
    "3S-hydroxyproline": r"O[C@@H]1[C@H](NCC1)C(=O)O",
    "3S-methyl-D-aspartic acid branched": r"C[C@@H]([C@H](C(=O)O)N)C(=O)O",
    "3S-methyl-D-aspartic acid wonky": r"N[C@H]([C@H](O)C(O)=O)C(O)=O",
    "3S-methylaspartic acid": r"C[C@@H]([C@@H](C(=O)O)N)C(=O)O",
    "3S-methylaspartic acid branched": r"C[C@@H]([C@@H](C(=O)O)N)C(=O)O",
    "3S-methylleucine": r"CC(C)[C@H](C)[C@@H](C(=O)O)N",
    "3S-methylproline": r"C[C@@H]1[C@H](NCC1)C(=O)O",
    "4,5-dehydroarginine": r"O=C(O)[C@@H](N)C/C=C/NC(N)=N",
    "4,5-dihydroxyornithine": r"C([C@@H](C(=O)O)N)C(C(N)O)O",
    "4-acetamidopyrrole-2-carboxylic acid": r"CC(=O)NC1=CNC(=C1)C(=O)O",
    "4-amino-2-hydroxy-3-isopropoxybenzoic acid": r"CC(C)OC1=C(C=CC(=C1O)C(=O)O)N",
    "4-aminobutyric acid": r"NCCCC(=O)O",
    "4-aminophenylalanine": r"C1=CC(=CC=C1C[C@@H](C(=O)O)N)N",
    "4-chlorobenzoic acid": r"C1=CC(=CC=C1C(=O)O)Cl",
    "4-hydroxy-3-nitrobenzoic acid": r"C1=CC(=C(C=C1C(=O)O)[N+](=O)[O-])O",
    "4-hydroxy-D-kynurenine": r"C1=C(O)C=C(C(=C1)C(=O)C[C@H](C(=O)O)N)N",
    "4-hydroxybenzoic acid": r"C1=CC(=CC=C1C(=O)O)O",
    "4-hydroxyglutamine": r"C(C(O)C(=O)N)[C@@H](C(=O)O)N",
    "4-hydroxyindole-3-carboxylic acid": r"c1cc2c(c(c1)O)c(c[nH]2)C(=O)O",
    "4-hydroxyphenylglycine": r"C1=CC(=CC=C1[C@@H](C(=O)O)N)O",
    "4-hydroxyphenylpyruvic acid": r"C1=CC(=CC=C1CC(=O)C(=O)O)O",
    "4-hydroxyproline": r"C1[C@H](NCC1O)C(=O)O",
    "4-hydroxythreonine": r"C([C@H]([C@@H](C(=O)O)N)O)O",
    "4-hydroxyvaline": r"CC(CO)[C@@H](C(=O)O)N",
    "4-methoxytryptophan": r"C1=CC=C2C(=C1OC)C(=CN2)C[C@@H](C(=O)O)N",
    "4-methylproline": r"CC1C[C@H](NC1)C(=O)O",
    "4-nitrotryptophan": r"C1=CC=C2C(=C1[N+](=O)[O-])C(=CN2)C[C@@H](C(=O)O)N",
    "4-oxoproline": r"C1[C@H](NCC1=O)C(=O)O",
    "4R-E-butenyl-4R-methylthreonine": r"C/C=C/C[C@@H](C)[C@H]([C@@H](C(=O)O)N)O",
    "4R-hydroxyproline": r"C1[C@H](NC[C@@H]1O)C(=O)O",
    "4R-methylproline": r"C[C@@H]1C[C@H](NC1)C(=O)O",
    "4R-propylproline": r"CCC[C@@H]1C[C@H](NC1)C(=O)O",
    "4S,5-dihydroxy-2S-aminopentanoic acid": r"O[C@@H](C[C@@H](C(=O)O)N)CO",
    "4S-acetyl-5S-methylproline": r"CC(=O)O[C@H]1C[C@H](N[C@H](C)1)C(=O)O",
    "4S-hydroxylysine": r"NCC[C@H](O)C[C@@H](C(=O)O)N",
    "4S-methylazetidine-2S-carboxylic acid": r"C[C@H]1C[C@H](N1)C(=O)O",
    "4S-methylproline": r"C[C@H]1C[C@H](NC1)C(=O)O",
    "4S-propenylproline": r"C/C=C\[C@H]1C[C@H](NC1)C(=O)O",
    "5,5-dimethylpipecolic acid": r"C1C(C)(C)CN[C@@H](C1)C(=O)O",
    "5-aminolevulinic acid": r"C(CC(=O)O)C(=O)CN",
    "5-chloroanthranilic acid": r"C1=CC(=C(C=C1Cl)C(=O)O)N",
    "5-chlorotryptophan": r"C1=CC2=C(C=C1Cl)C(=CN2)C[C@@H](C(=O)O)N",
    "5-methoxytyrosine": r"C1=C(OC)C(=CC=C1C[C@@H](C(=O)O)N)O",
    "5-methylorsellinic acid": r"C=1(C=C(C(=C(C1C)C)C(=O)O)O)O",
    "5S-methylproline": r"C1C[C@H](N[C@@H](C)1)C(=O)O",
    "6,7-dichlorotryptophan": r"C1=C(Cl)C(Cl)=C2C(=C1)C(=CN2)C[C@@H](C(=O)O)N",
    "6-chloro-4-hydroxy-1-methyl-indole-3-carboxylic acid": r"C(O)1=C(Cl)C=C2C(=C1)C(=CN(C)2)C(=O)O",
    "6-chloro-4-hydroxyindole-3-carboxylic acid": r"c(Cl)1cc2c(c(c1)O)c(c[nH]2)C(=O)O",
    "6-chlorotryptophan": r"C1=C(Cl)C=C2C(=C1)C(=CN2)C[C@@H](C(=O)O)N",
    "6-hydroxy-tetrahydro-isoquinoline-3-carboxylic acid": r"C1C(NCC2=C1C=C(C=C2)O)C(=O)O",
    "6-methylsalicylic acid": r"CC1=C(C(=CC=C1)O)C(=O)O",
    "6S-methyl-pipecolic acid": r"C1C[C@H](C)N[C@@H](C1)C(=O)O",
    "Acetyl-Coa": r"CC(=O)SCCNC(=O)CCNC(=O)C(C(C)(C)COP(=O)(O)OP(=O)(O)OC[C@@H]1[C@H]([C@H]([C@@H](O1)N2C=NC3=C(N=CN=C32)N)O)OP(=O)(O)O)O",
    "An acid hydrazine polyene (intermediate 14)": r"OC(=O)CCC(=O)NNCC(=O)O",
    "Compound 4 (formed by the decarboxylative condensation of L-Phe and succinyl-CoA)": r"C1=CC=C(C=C1)C[C@@H](C(=O)CCC(=O)O)N",
    "D-alanine": r"C[C@H](C(=O)O)N",
    "D-aspartic acid branched": r"C([C@H](C(=O)O)N)C(=O)O",
    "D-glutamic acid branched": r"C(CC(=O)O)[C@H](C(=O)O)N",
    "D-isovaline": r"CC[C@](C)(C(=O)O)N",
    "D-leucine": r"CC(C)C[C@H](C(=O)O)N",
    "D-lysergic acid": r"CN1C[C@@H](C=C2C1CC3=CNC4=CC=CC2=C34)C(=O)O",
    "D-phenylalanine": r"C1=CC=C(C=C1)C[C@H](C(=O)O)N",
    "D-phenyllactic acid": r"C1=CC=C(C=C1)C[C@H](C(=O)O)O",
    "D-pipecolic acid": r"C1CCN[C@H](C1)C(=O)O",
    "D-serine": r"C([C@H](C(=O)O)N)O",
    "Malonyl-CoA": r"CC(C)(COP(=O)(O)OP(=O)(O)OC[C@@H]1[C@H]([C@H]([C@@H](O1)N2C=NC3=C(N=CN=C32)N)O)OP(=O)(O)O)[C@H](C(=O)NCCC(=O)NCCSC(=O)CC(=O)O)O",
    "N-(1-methyl)-tryptophan": r"C1=CC=C2C(=C1)C(=CN(C)2)C[C@@H](C(=O)O)N",
    "N-(1-propargyl)-tryptophan": r"C1=CC=C2C(=C1)C(=CN(CC#C)2)C[C@@H](C(=O)O)N",
    "N-formylglycine": r"C(C(=O)O)NC=O",
    "N-hydroxyvaline": r"CC(C)[C@@H](C(=O)O)NO",
    "N-methylphenylalanine": r"CN[C@@H](CC1=CC=CC=C1)C(=O)O",
    "N-methyltyrosine": r"C1=CC(=CC=C1C[C@@H](C(=O)O)NC)O",
    "N1-methoxytryptophan": r"C1=CC=C2C(=C1)C(=CN(OC)2)C[C@@H](C(=O)O)N",
    "N5-acetyl-N5-hydroxyornithine": r"CC(=O)N(CCC[C@@H](C(=O)O)N)O",
    "N5-acetyl-hydroxyornithine": r"CC(=O)N(CCC[C@@H](C(=O)O)N)O",
    "N5-cis-anhydromevalonyl-N5-hydroxyornithine": r"C(N(C(=O)/C=C(/CCO)\C)O)CC[C@@H](C(O)=O)N",
    "N5-formyl-N5-hydroxyornithine": r"C(C[C@@H](C(=O)O)N)CN(C=O)O",
    "N5-hydroxyornithine": r"C(C[C@@H](C(=O)O)N)CNO",
    "N5-nitroso-N5-hydroxyornithine": r"O=NN(CCC[C@@H](C(=O)O)N)O",
    "N5-trans-anhydromevalonyl-N5-hydroxyornithine": r"C(C[C@@H](C(=O)O)N)CN(O)C(=O)/C=C(C)/CCO",
    "N6-hydroxylysine": r"C(CCNO)C[C@@H](C(=O)O)N",
    "O-dimethylallyl-L-tyrosine": r"CC(=CCOC1=CC=C(C=C1)C[C@@H](C(=O)O)N)C",
    "O-methylthreonine": r"C[C@H]([C@@H](C(=O)O)N)OC",
    "O-methyltyrosine": r"COC1=CC=C(C=C1)C[C@@H](C(=O)O)N",
    "R-3-hydroxy-3-methylproline": r"O[C@](C)1[C@H](NCC1)C(=O)O",
    "R-aza-beta-tyrosine": r"C1=CC(=NC=C1O)[C@@H](CC(=O)O)N",
    "R-beta-hydroxyphenylalanine": r"O[C@H](C1=CC=CC=C1)[C@@H](C(=O)O)N",
    "R-beta-hydroxytyrosine": r"C1=CC(=CC=C1[C@H]([C@@H](C(=O)O)N)O)O",
    "R-beta-methylphenylalanine": r"C[C@H](C1=CC=CC=C1)[C@@H](C(=O)O)N",
    "R-beta-methyltryptophan": r"C[C@H](C1=CNC2=CC=CC=C21)[C@@H](C(=O)O)N",
    "R-beta-phenylalanine": r"C1=CC=C(C=C1)[C@@H](CC(=O)O)N",
    "R-beta-tyrosine": r"C1=CC(=CC=C1[C@@H](CC(=O)O)N)O",
    "S-adenosylmethionine": r"C[S+](CC[C@@H](C(=O)[O-])N)C[C@@H]1[C@H]([C@H]([C@@H](O1)N2C=NC3=C(N=CN=C32)N)O)O",
    "S-beta-hydroxycyclohex-2S-enylalanine": r"C1C=C[C@H](CC1)[C@H](O)[C@@H](C(=O)O)N",
    "S-beta-hydroxyenduracididine": r"C1[C@H](NC(=N1)N)[C@H](O)[C@@H](C(=O)O)N",
    "S-beta-hydroxyphenylalanine": r"O[C@@H](C1=CC=CC=C1)[C@@H](C(=O)O)N",
    "S-beta-methylphenylalanine": r"C[C@@H](C1=CC=CC=C1)[C@@H](C(=O)O)N",
    "S-beta-tyrosine": r"C1=CC(=CC=C1[C@H](CC(=O)O)N)O",
    "acetic acid": r"CC(O)=O",
    "alanine": r"C[C@@H](C(=O)O)N",
    "alaninol": r"C[C@@H](CO)N",
    "allo-isoleucine": r"CC[C@@H](C)[C@@H](C(=O)O)N",
    "allo-threonine": r"C[C@@H]([C@@H](C(=O)O)N)O",
    "anthanillic acid": r"C1=CC=C(C(=C1)C(=O)O)N",
    "anthranilic acid": r"C1=CC=C(C(=C1)C(=O)O)N",
    "arginine": r"C(C[C@@H](C(=O)O)N)CN=C(N)N",
    "argininol": r"N[C@H](CO)CCCN=C(N)N",
    "asparagine": r"C([C@@H](C(=O)O)N)C(=O)N",
    "aspartic acid": r"C([C@@H](C(=O)O)N)C(=O)O",
    "aspartic acid branched": r"C([C@@H](C(=O)O)N)C(=O)O",
    "azetidine-2-carboxylic acid": r"O=C(O)[C@H]1NCC1",
    "benzoic acid": r"C1=CC=C(C=C1)C(=O)O",
    "benzoxazolinate": r"c1ccc2c(c1)nc(o2)C(=O)O",
    "beta-alanine": r"NCCC(=O)O",
    "beta-hydroxy-3-hydroxy-O-methyl-5-methyltyrosine": r"C1=C(C)C(=C(O)C=C1C(O)[C@@H](C(=O)O)N)OC",
    "beta-hydroxy-gamma-methyl-hexadecanoic acid": r"CCCCCCCCCCCCC(C)C(O)CC(=O)O",
    "beta-hydroxyarginine": r"C(C(O)[C@@H](C(=O)O)N)CN=C(N)N",
    "beta-hydroxyphenylalanine": r"OC(C1=CC=CC=C1)[C@@H](C(=O)O)N",
    "beta-hydroxytyrosine": r"C1=CC(=CC=C1C([C@@H](C(=O)O)N)O)O",
    "beta-lysine": r"C(C[C@@H](CC(=O)O)N)CN",
    "beta-methylphenylalanine": r"CC(C1=CC=CC=C1)C(C(=O)O)N",
    "beta-tyrosine": r"C1=CC(=CC=C1C(CC(=O)O)N)O",
    "betaine": r"C[N+](C)(C)CC(=O)O",
    "butyric acid": r"CCCC(=O)O",
    "caffeic acid": r"OC(=O)\C=C\c1ccc(O)c(O)c1",
    "capreomycidine": r"C1CN=C(N[C@H]1[C@@H](C(=O)O)N)N",
    "cinnamic acid": r"C1=CC=C(C=C1)/C=C/C(=O)O",
    "citrulline": r"C(C[C@@H](C(=O)O)N)CNC(=O)N",
    "colletorin D acid": r"CC1=CC(=C(C(=C1C(=O)O)O)CC=C(C)C)O",
    "coumaric acid": r"C1=CC(=CC=C1/C=C/C(=O)O)O",
    "cysteic acid": r"C([C@@H](C(=O)O)N)S(=O)(=O)O",
    "cysteine": r"C([C@@H](C(=O)O)N)S",
    "cysteine branched": r"C([C@@H](C(=O)O)N)S",
    "dehydroarginine": r"C(CN=C(N)N)/C=C(/C(=O)O)\N",
    "dehydrophenylalanine": r"N/C(=C\C1=CC=CC=C1)/C(=O)O",
    "dehydrotryptophan": r"C1=CC=C2C(=C1)C(=CN2)/C=C(/C(=O)O)\N",
    "dehydrovaline": r"CC(=C(C(=O)O)N)C",
    "dihydrolysergic acid": r"CN1CC(CC2C1CC3=CNC4=CC=CC2=C34)C(=O)O",
    "dimethylsulfoniopropionic acid": r"C[S+](C)CCC(=O)O",
    "enduracididine": r"C1[C@H](NC(=N1)N)C[C@@H](C(=O)O)N",
    "fatty acid": r"O=C(O)C*",
    "fumaric acid": r"C(=C/C(=O)O)\C(=O)O",
    "glutamic acid": r"C(CC(=O)O)[C@@H](C(=O)O)N",
    "glutamine": r"C(CC(=O)N)[C@@H](C(=O)O)N",
    "glycine": r"NCC(=O)O",
    "glycolic acid": r"C(C(=O)O)O",
    "graminine": r"O=NN(O)CCC[C@H](N)(C(=O)O)",
    "grifolic acid": r"CC(C)=CCC/C(C)=C/CC/C(C)=C/CC1=C(O)C=C(C)C(C(=O)O)=C(O)1",
    "guanidinoacetic acid": r"C(C(=O)O)N=C(N)N",
    "histidine": r"C1=C(NC=N1)C[C@@H](C(=O)O)N",
    "homophenylalanine": r"C1=CC=C(C=C1)CC[C@@H](C(=O)O)N",
    "homoserine": r"C(CO)[C@@H](C(=O)O)N",
    "homotyrosine": r"C1=CC(=CC=C1CC[C@@H](C(=O)O)N)O",
    "hydroxyproline": r"C(*)1C[C@H](NC(*)1)C(=O)O",
    "indole pyruvic acid": r"C1=CC=C2C(=C1)C(=CN2)CC(=O)C(=O)O",
    "isoleucine": r"CC[C@H](C)[C@@H](C(=O)O)N",
    "isovaline": r"CC[C@@](C)(C(=O)O)N",
    "kynurenine": r"C1=CC=C(C(=C1)C(=O)C[C@@H](C(=O)O)N)N",
    "lactic acid": r"CC(C(=O)O)O",
    "leucine": r"CC(C)C[C@@H](C(=O)O)N",
    "leucinol": r"CC(C)C[C@@H](CO)N",
    "linoleic acid": r"CCCCC/C=C\C/C=C\CCCCCCCC(=O)O",
    "lysine": r"C(CCN)C[C@@H](C(=O)O)N",
    "malic acid": r"C(C(C(=O)O)O)C(=O)O",
    "malonamate": r"NC(=O)CC(=O)O",
    "meta-tyrosine": r"C1=CC(=CC(=C1)O)C[C@@H](C(=O)O)N",
    "methionine": r"CSCC[C@@H](C(=O)O)N",
    "methylglutaconyl hydroxyornithine": r"C/C(=C\C(=O)N(CCCC(C(=O)O)N)O)/CC(=O)O",
    "nicotinic acid": r"C1=CC(=CN=C1)C(=O)O",
    "norcoronamic acid": r"C[C@H]1C[C@]1(C(=O)O)N",
    "ochratoxin beta": r"CC1CC2=C(C(=C(C=C2)C(=O)O)O)C(=O)O1",
    "ornithine": r"C(C[C@@H](C(=O)O)N)CN",
    "p-hydroxybenzoylformic acid": r"C1=CC(=CC=C1C(=O)C(=O)O)O",
    "p-hydroxymandelate": r"C1=CC(=CC=C1C(C(=O)O)O)O",
    "para-aminobenzoic acid": r"O=C(O)c1ccc(N)cc1",
    "pentanoic acid": r"CCCCC(=O)O",
    "phenazine-1,6-dicarboxylic acid": r"C1=CC(=C2C(=C1)N=C3C(=N2)C=CC=C3C(=O)O)C(=O)O",
    "phenylalanine": r"C1=CC=C(C=C1)C[C@@H](C(=O)O)N",
    "phenylalaninol": r"C1=CC=C(C=C1)C[C@@H](CO)N",
    "phenylglycine": r"C1=CC=C(C=C1)[C@@H](C(=O)O)N",
    "phenyllactic acid": r"C1=CC=C(C=C1)C[C@@H](C(=O)O)O",
    "phenylpyruvic acid": r"C1=CC=C(C=C1)CC(=O)C(=O)O",
    "pipecolic acid": r"C1CCN[C@@H](C1)C(=O)O",
    "piperazic acid": r"C1C[C@H](NNC1)C(=O)O",
    "piperonylic acid": r"OC(=O)c1ccc2OCOc2c1",
    "proline": r"C1C[C@H](NC1)C(=O)O",
    "propionic acid": r"CCC(=O)O",
    "pyrrole-2-carboxylic acid": r"C1=CNC(=C1)C(=O)O",
    "pyruvic acid": r"CC(=O)C(=O)O",
    "quinoxaline-2-carboxylic acid": r"C1=CC=C2C(=C1)N=CC(=N2)C(=O)O",
    "salicylic acid": r"C1=CC=C(C(=C1)C(=O)O)O",
    "serine": r"C([C@@H](C(=O)O)N)O",
    "succinic semialdehyde": r"C(CC(=O)O)C=O",
    "succinyl-hydrazinoacetic acid": r"N/N=C/C=C/C=C/C=C/C=C/C=C/C(=O)O",
    "tetradecanoic acid": r"CCCCCCCCCCCCCC(=O)O",
    "threonine": r"C[C@H]([C@@H](C(=O)O)N)O",
    "trans-2-crotylglycine": r"C/C=C/C[C@@H](C(=O)O)N",
    "trans-2-hexenoic acid": r"CCC\C=C\C(O)=O",
    "tricarballylic acid": r"C(C(CC(=O)O)C(=O)O)C(=O)O",
    "tryptophan": r"C1=CC=C2C(=C1)C(=CN2)C[C@@H](C(=O)O)N",
    "tyrosine": r"C1=CC(=CC=C1C[C@@H](C(=O)O)N)O",
    "ustethylinic acid": r"c1(C)c(O)c(C(=O)O)c(CC)cc(O)1",
    "valine": r"CC(C)[C@@H](C(=O)O)N",
    "valine isocyanide": r"CC(C)[C@H]([N+]#[C-])C(O)=O",
    "valinol": r"CC(C)[C@@H](CO)N",
}


@dataclass
class ADomain:
    """
    Dataclass representing an A domain.
    
    :param protein: Name of the protein containing the A domain.
    :param start: Start position of the A domain.
    :param end: End position of the A domain.
    :param domain_nr: Domain number of A domain in NRPS (optional).
    :param sequence: Amino acid sequence of the A domain (optional).
    :param extended_signature: Extended signature of the A domain (optional).
    """

    protein: str
    start: int
    end: int
    domain_nr: int | None = None
    sequence: str | None = None
    extended_signature: str | None = None


def _b2s(x: Any) -> str:
    """
    Convert input to string.
    
    :param x: input object
    :return: string representation
    """
    if isinstance(x, (bytes, bytearray)):
        return x.decode()

    if hasattr(x, "sequence"):
        s = x.sequence
        return s.decode() if isinstance(s, (bytes, bytearray)) else str(s)
    
    return str(x)


def extract_domain_hits(
    seq_id: str,
    sequence: str,
    evalue_cutoff: float = 1e-5,
) -> list[dict[str, Any]]:
    """
    Extract domain hits from a given protein sequence using HMMER.

    :param seq_id: Identifier for the protein sequence.
    :param sequence: Amino acid sequence of the protein.
    :param evalue_cutoff: E-value cutoff for HMMER hits.
    :return: List of dictionaries representing domain hits.
    :raises ImportError: If HMMER is not available.
    """
    if not PYHMMER_AVAILABLE:
        raise ImportError("HMMER is not available. Please install pyhmmer to use this function.")

    alphabet = easel.Alphabet.amino()
    text_seq = easel.TextSequence(name=seq_id.encode(), sequence=sequence)
    seq = text_seq.digitize(alphabet)

    hits_iter = hmmer.hmmscan([seq], HMM_DB, cpus=1, E=evalue_cutoff)

    query_hits = next(hits_iter)  # expect only one sequence

    out = []
    for hit in query_hits:
        model_name = _b2s(hit.name)

        for dom in hit.domains:
            q_from = int(dom.env_from)
            q_to = int(dom.env_to)

            aln = dom.alignment
            hmm_aln = _b2s(aln.hmm_sequence)
            query_aln = _b2s(aln.target_sequence)

            out.append(
                dict(
                    seq_id=seq_id,
                    model=model_name,
                    q_from=q_from,
                    q_to=q_to,
                    evalue=float(dom.i_evalue),
                    score=float(dom.score),
                    hmm_aln=hmm_aln,
                    query_aln=query_aln,
                    domain_obj=dom,
                )
            )

    out.sort(key=lambda d: (d["q_from"], d["q_to"], d["model"]))

    return out


def pair_domains(
    domain_hits: list[dict[str, Any]],
    max_gap: int = 200,
) -> list[tuple[ADomain, str, str]]:
    """
    Pair AMP-binding and AMP-binding_C domain hits.
    
    :param domain_hits: List of domain hit dictionaries.
    :param max_gap: Maximum allowed gap between paired domains.
    :return: List of tuples containing ADomain objects and their alignments.
    """
    hits = sorted(domain_hits, key=lambda d: d["q_from"])

    a_domains: list[ADomain] = []
    for h1 in hits:
        if h1["model"] != "AMP-binding":
            continue
        
        n_from, n_to = h1["q_from"], h1["q_to"]

        matched = None
        for h2 in hits:
            if h2["model"] != "AMP-binding_C":
                continue

            c_from = h2["q_from"]

            if c_from > n_to and (c_from - n_to) <= max_gap:
                matched = h2
                break

        start0 = n_from - 1
        end0 = matched["q_to"] if matched is not None else n_to
        a_domains.append((ADomain(
            protein=h1["seq_id"],
            start=start0,
            end=end0),
            h1["hmm_aln"],
            h1["query_aln"]
        ))

    a_domains.sort(key=lambda t: t[0].start)
    for i, (d, _, _) in enumerate(a_domains, start=1):
        d.domain_nr = i
    
    return a_domains


def extract_signature_from_alignment(hmm_aln: str, query_aln: str) -> str | None:
    """
    Extract the extended signature from the given HMM and query alignments.
    
    :param hmm_aln: HMM alignment string
    :param query_aln: query alignment string
    :return: extended signature string or None if invalid
    """
    wanted = set(POSITIONS_ACTIVE_SITE)
    picked: dict[int, str] = {}

    hmm_pos = 0  # 1-based counter, increment when HMM char is not a gap

    for h, q in zip(hmm_aln, query_aln):
        if h != "-":
            hmm_pos += 1
            if hmm_pos in wanted and hmm_pos not in picked:
                picked[hmm_pos] = q

    # Quick fix
    missing = wanted - set(picked.keys())
    for m in missing:
        picked[m] = "-"

    out = []
    for p in POSITIONS_ACTIVE_SITE:
        if p not in picked:
            return None
        out.append(picked[p])
    
    sig = "".join(out).upper()
    if not sig or not all(c in VALID for c in sig):
        return None
    
    return sig


def fill_domain_sequences(
    domains: list[ADomain],
    protein_seq: str,
    min_len: int = 100,
) -> list[ADomain]:
    """
    Fill in the sequences for the given domains from the protein sequence.

    :param domains: List of ADomain objects.
    :param protein_seq: Amino acid sequence of the protein.
    :param min_len: Minimum length of domain sequence to keep.
    :return: List of ADomain objects with sequences filled in.
    """
    out = []

    for d in domains:
        seq = protein_seq[d.start:d.end]
        if len(seq) >= min_len:
            d.sequence = seq
            out.append(d)

    return out


def find_a_domains(
    seq_id: str,
    protein_seq: str,
    evalue_cutoff: float = 1e-5,
) -> list[ADomain]:
    """
    Find A domains in a given protein sequence using HMMER.
    
    :param seq_id: Identifier for the protein sequence.
    :param protein_seq: Amino acid sequence of the protein.
    :param evalue_cutoff: E-value cutoff for HMMER hits.
    :return: List of ADomain objects representing found A domains.
    """
    hits = extract_domain_hits(seq_id, protein_seq, evalue_cutoff)

    hits = [h for h in hits if h["model"] in {"AMP-binding", "AMP-binding_C"}]

    paired = pair_domains(hits, max_gap=200)

    domains_only: list[ADomain] = []
    for d, hmm_aln, query_aln in paired:
        d.extended_signature = extract_signature_from_alignment(hmm_aln, query_aln)
        domains_only.append(d)

    domains_only = fill_domain_sequences(domains_only, protein_seq, min_len=100)

    domains_only = [d for d in domains_only if d.extended_signature is not None]

    domains_only.sort(key=lambda d: (d.protein, d.start))

    return domains_only


def featurize_signature(sig: str) -> np.ndarray:
    """
    Featurize the given extended signature into a numerical feature array.

    :param sig: Extended signature string.
    :return: Numpy array of features.
    """
    assert len(sig) == len(POSITIONS_ACTIVE_SITE), "signature length mismatch"

    features: np.ndarray = np.zeros((len(POSITIONS_ACTIVE_SITE), len(FEATURE_NAMES)), dtype=np.float32)
    for i, aa in enumerate(sig):
        aa_feats = FEATURES.get(aa)
        if aa_feats is None:
            raise ValueError(f"invalid amino acid '{aa}' in signature")
        features[i, :] = np.array(aa_feats, dtype=np.float32)
    
    return features.flatten()  # shape (n_positions * n_features,)


def get_paras_model_path(cache_dir: Path) -> Path:
    """
    Get the path to the cached PARAS model, downloading it if necessary.

    :param cache_dir: Path to the cache directory.
    :return: Path to the PARAS model file.
    """
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)

    # If we already know the path and it still exists, reuse it
    cached = _PARAS_MODEL_PATH_CACHE.get(PARAS_DOWNLOAD_URL)
    if cached is not None and cached.exists():
        return cached
    
    # Otherwise download/prepare and remember the path
    model_path = download_and_prepare(PARAS_DOWNLOAD_URL, cache_dir)
    _PARAS_MODEL_PATH_CACHE[PARAS_DOWNLOAD_URL] = model_path
    return model_path


def resolve_paras_model_path(model_path: Path | None, cache_dir: Path) -> Path:
    """
    Resolve the PARAS model path, downloading it to cache if necessary.

    :param model_path: User-specified model path or None.
    :param cache_dir: Path to the cache directory.
    :return: Path to the PARAS model file.
    """
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)

    if model_path is not None:
        model_path = Path(model_path).expanduser().resolve()
        if not model_path.exists():
            raise FileNotFoundError(f"specified PARAS model path does not exist: {model_path}")
        return model_path
    
    # No model given -> download/prepare in cache dir
    return Path(download_and_prepare(PARAS_DOWNLOAD_URL, cache_dir))


def load_paras_model(model_path: Path) -> object:
    """
    Load the PARAS model from the given path.

    :param model_path: Path to the PARAS model file.
    :return: Loaded RandomForestClassifier model.
    """
    key = str(model_path)
    model = _PARAS_MODEL_OBJ_CACHE.get(key)
    
    if model is None:
        log.info(f"loading PARAS model: {model_path}")
        model = joblib.load(model_path)
        _PARAS_MODEL_OBJ_CACHE[key] = model
    
    return model


class ParasModel(DomainInferenceModel):
    """
    Model for predicting A domain substrate specificity using PARAS.

    :param domain: The domain to make predictions for.
    :param threshold: Probability threshold for predictions.
    :return: A list of InferenceResult objects containing the predictions.
    """

    name: str = "paras"

    def __init__(
        self,
        threshold: float = 0.1,
        keep_top: int = 3,
        cache_dir: Path | str | None = None,
        model_path: Path | str | None = None,
    ) -> None:
        """
        Initialize the ParasModel.

        :param threshold: Probability threshold for predictions.
        :param keep_top: Number of top predictions to keep.
        :param cache_dir: Directory to cache the model.
        :param model_path: Path to a custom PARAS model file.
        """
        super().__init__()
        
        # Set cache directory
        if cache_dir is None:
            cache_dir = PARAS_CACHE_DIR
        self.cache_dir = Path(cache_dir)
        self.model_path = Path(model_path) if model_path is not None else None

        # Set other parameters
        self.threshold = threshold
        self.keep_top = keep_top

    def __post_init__(self) -> None:
        """
        Post-initialization to set up the model.

        :raises ValueError: If threshold is not between 0 and 1.
        """
        super().__post_init__()

        # Make sure threshold is between 0 and 1
        if not (0.0 <= self.threshold <= 1.0):
            raise ValueError("threshold must be between 0 and 1")
        
        # Make sure keep_top is an int and at least 1
        if not (isinstance(self.keep_top, int) and self.keep_top >= 1):
            raise ValueError("keep_top must be an integer >= 1")


    def predict(self, domain: Domain) -> list[InferenceResult]:
        """
        Predict the substrate specificity for the given domain.
        
        :param domain: The domain to make predictions for.
        :return: A list of InferenceResult objects containing the predictions.
        """
        if domain.type == "AMP-binding":
            # Prepare model
            model_file = resolve_paras_model_path(self.model_path, self.cache_dir)
            model = load_paras_model(model_file)

            # Find A domains in the sequence
            a_domains = find_a_domains(seq_id=domain.id, protein_seq=domain.sequence)

            # Make predictions
            unknown_prediction = self.result(label="unknown", score=0.0, metadata={})

            match a_domains:
                case []:
                    log.warning(f"no A domains found in sequence {domain.id}; unable to predict substrate")
                    return [unknown_prediction]
                case [a_domain]:
                    sig = a_domain.extended_signature
                    features = featurize_signature(sig).reshape(1, -1)
                    pred = model.predict_proba(features)
                    
                    # Identify top predictions
                    lbls = model.classes_
                    top_indices = np.argsort(pred, axis=1)[0][-self.keep_top:][::-1]
                    top_lbls = [lbls[i] for i in top_indices]
                    top_prbs = [pred[0, i] for i in top_indices]

                    results = []
                    for top_lbl, top_prb in zip(top_lbls, top_prbs):
                        if top_prb >= self.threshold:
                            smi = LABEL_TO_SMILES.get(top_lbl)
                            metadata = {}
                            if smi is None:
                                log.warning(f"no SMILES found for predicted label '{top_lbl}'; returning label only")
                            else:
                                metadata["smiles"] = smi
                            results.append(self.result(label=top_lbl, score=round(float(top_prb), 4), metadata=metadata))
                        else:
                            log.debug(f"prediction '{top_lbl}' for sequence {domain.id} below threshold ({top_prb:.4f} < {self.threshold}); skipping")

                    if not results:
                        return [unknown_prediction]

                    return results
                case _:
                    log.error(f"found multiple ({len(a_domains)}) A domains in sequence {domain.id}; unable to predict substrate")
                    return [unknown_prediction]

        else:
            # Not the domain type of interest
            return []
