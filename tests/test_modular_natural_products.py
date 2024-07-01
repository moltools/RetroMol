"""Integration tests for parsing and sequencing modular natural products."""

import json
import os
import typing as ty

import pytest

from retromol.retrosynthesis.chem import Molecule
from retromol.retrosynthesis.helpers import RetroMolTimeoutError, timeout
from retromol.retrosynthesis.parsing import (
    Result,
    parse_mol,
    parse_molecular_patterns,
    parse_reaction_rules,
)

TESTS_DIR_NAME = os.path.dirname(__file__)
# PATH_TO_MONOMERS = os.path.join(TESTS_DIR_NAME, "fixtures", "monomers.json")
PATH_TO_MONOMERS = os.path.join(os.path.dirname(TESTS_DIR_NAME), "app", "src", "server", "data", "monomers.json")
# PATH_TO_REACTIONS = os.path.join(TESTS_DIR_NAME, "fixtures", "reactions.json")
PATH_TO_REACTIONS = os.path.join(os.path.dirname(TESTS_DIR_NAME), "app", "src", "server", "data", "reactions.json")
MONOMERS_SRC = json.load(open(PATH_TO_MONOMERS, "r", encoding="utf-8"))
REACTIONS_SRC = json.load(open(PATH_TO_REACTIONS, "r", encoding="utf-8"))
MONOMERS = parse_molecular_patterns(MONOMERS_SRC)
REACTIONS = parse_reaction_rules(REACTIONS_SRC)


@timeout(10)
def parse_mol_timed(mol: Molecule) -> Result:
    """Parse a molecule with a timeout.

    :param mol: The molecule to parse.
    :type mol: Molecule
    :return: The result of the parsing.
    :rtype: Result
    """
    try:
        result = parse_mol(mol, REACTIONS, MONOMERS)
    except Exception:
        result = Result(mol.name, mol.compiled, False)

    return result


@pytest.mark.parametrize(
    "compound_name, test, expected, expected_raises",
    [
        (
            "10-deoxymethynolide",
            r"CC[C@@H]1[C@@H](/C=C/C(=O)[C@@H](C[C@@H]([C@@H]([C@H](C(=O)O1)C)O)C)C)C",
            [
                "polyketide|B2",
                "polyketide|C1",
                "polyketide|A2",
                "polyketide|D2",
                "polyketide|B2",
            ],
            None,
        ),
        (
            "13-deoxytedanolide",
            r"C/C=C\[C@H](C)[C@@H]1[C@](O1)(C)[C@H]([C@H]2COC(=O)[C@@H]([C@H]([C@@H](C(=O)[C@@H]([C@H](/C(=C/[C@@H](C(=O)CC[C@H](C2=O)C)C)/C)O)C)C)OC)O)O",
            [
                "polyketide|C2",
                "polyketide|C2",
                "polyketide|B6",
                "polyketide|A2",
                "polyketide|D1",
                "polyketide|A2",
                "polyketide|C2",
                "polyketide|B2",
                "polyketide|A2",
                "polyketide|B5",
            ],
            None,
        ),
        (
            "6-deoxyerythronolide_B",
            r"CC[C@@H]1[C@@H]([C@@H]([C@H](C(=O)[C@@H](C[C@@H]([C@@H]([C@H]([C@@H]([C@H](C(=O)O1)C)O)C)O)C)C)C)O)C",
            [
                "polyketide|B2",
                "polyketide|B2",
                "polyketide|A2",
                "polyketide|D2",
                "polyketide|B2",
                "polyketide|B2",
            ],
            None,
        ),
        (
            "abyssomicin_C",
            r"C[C@@H]1C[C@]23OC(=O)C4=C2OC1[C@H](O)C3\C=C\C(=O)[C@@H](C)C[C@@H](C)C4=O",
            [],
            None,
        ),
        (
            "acutiphycin",
            r"CCCCC[C@@H]1C/C=C(\[C@H](C(C(=O)[C@H](/C=C(\[C@@H]2C[C@@H](C[C@@](O2)(CC(=O)O1)O)O)/C)C)(C)C)O)/C",
            [
                "polyketide|D1",
                "polyketide|D1",
                "polyketide|B1",
                "polyketide|C2",
                "polyketide|B3",
                "polyketide|A2",
                "polyketide|C2",
                "polyketide|B1",
                "polyketide|B1",
                "polyketide|A1",
            ],
            None,
        ),
        (
            "amicoumacin",
            r"CC(C)C[C@@H]([C@@H]1CC2=C(C(=CC=C2)O)C(=O)O1)NC(=O)[C@H]([C@H]([C@H](CC(=O)N)N)O)O",
            [
                "peptide|pubchem|236",
                "polyketide|B5",
                "peptide|pubchem|857",
                "polyketide|B1",
                "polyketide|C1",
                "polyketide|C1",
                "polyketide|A1",
            ],
            None,
        ),
        (
            "amphidinolide_J",
            r"CCC/C=C/[C@@H](C)[C@H]1C(/C=C\C([C@H](C=CCCC(=C)[C@H](CC(=O)O1)C)O)C)O",
            [
                "polyketide|D1",
                "polyketide|C1",
                "polyketide|D9",
                "polyketide|B1",
                "polyketide|C2",
                "polyketide|B1",
                "polyketide|C1",
                "polyketide|D12",
                "polyketide|D8",
            ],
            None,
        ),
        (
            "amphidinolide_P",
            r"C[C@@H]1C(=C)C[C@H]2[C@H]3[C@@H](O3)CC(=C)/C=C/[C@H](OC(=O)C[C@@]1(O2)O)[C@H](C)C(=C)C",
            [
                "polyketide|D9",
                "polyketide|C1",
                "polyketide|A8",
                "polyketide|C1",
                "polyketide|B1",
                "polyketide|A9",
                "polyketide|A1",
            ],
            None,
        ),
        (
            "anthracimycin",
            r"C[C@@H]1/C=C\C=C\[C@H](OC(=O)[C@@H](C(=O)/C=C(/[C@H]2[C@@H]1C=C[C@@H]3[C@@H]2CC=C(C3)C)\O)C)C",
            [
                "polyketide|B1",
                "polyketide|C1",
                "polyketide|C2",
                "polyketide|C1",
                "polyketide|C1",
                "polyketide|D2",
                "polyketide|C1",
                "polyketide|C1",
                "polyketide|A1",
                "polyketide|A2",
            ],
            None,
        ),
        (
            "bitungolide_F",
            r"CC[C@@H]1C=CC(=O)O[C@@H]1[C@H](C)CC[C@H](C[C@@H](/C=C/C=C/C2=CC=CC=C2)O)O",
            [
                "polyketide|C1",
                "polyketide|C1",
                "polyketide|B1",
                "polyketide|B1",
                "polyketide|D2",
                "polyketide|B4",
                "polyketide|C1",
            ],
            None,
        ),
        (
            "borrelidin",
            r"C[C@H]1C[C@H](C[C@@H]([C@H](/C(=C\C=C\C[C@H](OC(=O)C[C@@H]([C@H](C1)C)O)[C@@H]2CCC[C@H]2C(=O)O)/C#N)O)C)C",
            [
                "polyketide|B1",
                "polyketide|C1",
                "polyketide|C13",
                "polyketide|B2",
                "polyketide|D2",
                "polyketide|D2",
                "polyketide|D2",
                "polyketide|B1",
            ],
            None,
        ),
        (
            "callystatin_A",
            r"CC[C@H](C)[C@H]([C@H](C)C(=O)[C@H](C)/C=C(\C)/C=C/C[C@@H](C)/C=C(/CC)\C=C\[C@H]1CC=CC(=O)O1)O",
            [
                "polyketide|D2",
                "polyketide|B2",
                "polyketide|A2",
                "polyketide|C2",
                "polyketide|C1",
                "polyketide|D2",
                "polyketide|C4",
                "polyketide|C1",
                "polyketide|B1",
                "polyketide|C1",
            ],
            None,
        ),
        # (
        #     "carolacton",
        #     r"C[C@@H]\1CCC[C@@H]([C@H](OC(=O)[C@@H]([C@@H](/C=C1)O)O)/C(=C/[C@@H](C)C(=O)[C@H](C)[C@@H](CC(=O)O)OC)/C)C",
        #     [],
        #     None,
        # ),
        (
            "chaetoglobosin_A",
            r"C[C@H]\1C/C=C/[C@H]2[C@H]3[C@](O3)([C@H]([C@@H]4[C@@]2(C(=O)/C=C/C(=O)[C@@H](/C(=C1)/C)O)C(=O)N[C@H]4CC5=CNC6=CC=CC=C65)C)C",
            [],
            None,
        ),
        (
            "chlorotonil_A",
            r"C[C@@H]1/C=C\C=C\[C@@H](OC(=O)[C@H](C(=O)C(C(=O)[C@@H]2[C@H]1C=C[C@H]3[C@H]2[C@@H](C=C(C3)C)C)(Cl)Cl)C)C",
            [
                "polyketide|B1",
                "polyketide|C1",
                "polyketide|C2",
                "polyketide|C1",
                "polyketide|C1",
                "polyketide|D2",
                "polyketide|C2",
                "polyketide|C1",
                "polyketide|A1",
                "polyketide|A2",
            ],
            None,
        ),
        (
            "daptomycin",
            r"CCCCCCCCCC(=O)NC(CC1=CNC2=CC=CC=C21)C(=O)NC(CC(=O)N)C(=O)NC(CC(=O)O)C(=O)NC3C(OC(=O)C(NC(=O)C(NC(=O)C(NC(=O)CNC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)CNC3=O)CCCN)CC(=O)O)C)CC(=O)O)CO)C(C)CC(=O)O)CC(=O)C4=CC=CC=C4N)C",
            [
                "polyketide|D1",
                "polyketide|D1",
                "polyketide|D1",
                "polyketide|D1",
                "peptide|pubchem|1148",
                "peptide|pubchem|236",
                "peptide|pubchem|424",
                "peptide|pubchem|205",
                "peptide|pubchem|750",
                "peptide|pubchem|389",
                "peptide|pubchem|424",
                "peptide|pubchem|602",
                "peptide|pubchem|424",
                "peptide|pubchem|750",
                "peptide|pubchem|617",
                "peptide|pubchem|237657",
                "peptide|pubchem|846",
            ],
            None,
        ),
        (
            "dictyostatin",
            r"C[C@H]1CC[C@H]([C@@H]([C@@H](OC(=O)/C=C\C=C\[C@H]([C@H](C[C@@H](/C=C\[C@@H]([C@@H]([C@H](C1)C)O)C)O)O)C)[C@@H](C)/C=C\C=C)C)O",
            [
                "polyketide|C1",
                "polyketide|C2",
                "polyketide|B2",
                "polyketide|B1",
                "polyketide|D2",
                "polyketide|D2",
                "polyketide|B2",
                "polyketide|C1",
                "polyketide|B1",
                "polyketide|B2",
                "polyketide|C1",
                "polyketide|C1",
            ],
            None,
        ),
        (
            "discodermolide",
            r"C[C@H]1[C@@H](OC(=O)[C@@H]([C@H]1O)C)C[C@@H](/C=C\[C@H](C)[C@@H]([C@@H](C)/C=C(/C)\C[C@H](C)[C@H]([C@H](C)[C@H]([C@@H](C)/C=C\C=C)OC(=O)N)O)O)O",
            [
                "polyketide|C1",
                "polyketide|C2",
                "polyketide|B2",
                "polyketide|B2",
                "polyketide|D2",
                "polyketide|C2",
                "polyketide|B2",
                "polyketide|C1",
                "polyketide|B1",
                "polyketide|B2",
                "polyketide|B2",
            ],
            None,
        ),
        (
            "epothilone",
            r"C[C@H]1CCC[C@@H]2[C@@H](O2)C[C@H](OC(=O)C[C@H](C(C(=O)[C@@H]([C@H]1O)C)(C)C)O)/C(=C/C3=CSC(=N3)C)/C",
            [
                "peptide|pubchem|594",
                "polyketide|C2",
                "polyketide|B1",
                "polyketide|C1",
                "polyketide|D1",
                "polyketide|D2",
                "polyketide|B2",
                "polyketide|A3",
                "polyketide|B1",
            ],
            None,
        ),
        (
            "erythromycin",
            r"CCC1C(C(C(C(=O)C(CC(C(C(C(C(C(=O)O1)C)OC2CC(C(C(O2)C)O)(C)OC)C)OC3C(C(CC(O3)C)N(C)C)O)(C)O)C)C)O)(C)O",
            [
                "polyketide|B7",
                "polyketide|B2",
                "polyketide|A2",
                "polyketide|D7",
                "polyketide|B2",
                "polyketide|B2",
            ],
            None,
        ),
        (
            "georatusin",
            r"CC[C@@H](C)C[C@H](C)[C@@H]1[C@H](C[C@@H]([C@@H]2[C@H](C[C@@H]([C@@](O2)([C@H](C(=O)N[C@@H](C(=O)O1)CC3=CNC4=CC=CC=C43)C)O)C)C)C)C",
            [
                "polyketide|D2",
                "polyketide|D2",
                "polyketide|B2",
                "polyketide|D2",
                "polyketide|B2",
                "polyketide|D2",
                "polyketide|A2",
                "peptide|pubchem|1148",
            ],
            None,
        ),
        (
            "gephyronic_acid",
            r"C[C@@H]1[C@@H](O[C@@](C([C@H]1OC)(C)C)([C@@H](C)C[C@H](C)[C@@H]([C@@]2([C@H](O2)[C@@H](C)C=C(C)C)C)O)O)CC(=O)O",
            [
                "polyketide|C2",
                "polyketide|C2",
                "polyketide|B2",
                "polyketide|D2",
                "polyketide|A3",
                "polyketide|B2",
                "polyketide|B1",
            ],
            None,
        ),
        (
            "harzianic_acid",
            r"CCC/C=C/C=C/C(=C\1/C(=O)C(N(C1=O)C)CC(C(C)C)(C(=O)O)O)/O",
            [
                "polyketide|D1",
                "polyketide|C1",
                "polyketide|C1",
                "polyketide|A1",
                "peptide|pubchem|29",
            ],
            None,
        ),
        (
            "herboxidiene",
            r"C[C@H]1CC[C@@H](O[C@@H]1/C(=C/C=C/[C@@H](C)C[C@@]2([C@H](O2)[C@H](C)[C@H]([C@@H](C)O)OC)C)/C)CC(=O)O",
            [
                "polyketide|B2",
                "polyketide|C2",
                "polyketide|D2",
                "polyketide|C1",
                "polyketide|C2",
                "polyketide|B2",
                "polyketide|D1",
                "polyketide|C1",
            ],
            None,
        ),
        (
            "hymenosetin",
            r"C/C=C/[C@@H]1C(=C[C@@H]2C[C@@H](CC[C@H]2[C@]1(C)/C(=C\3/C(=O)[C@@H](NC3=O)[C@@H](C)O)/O)C)C",
            [
                "polyketide|C1",
                "polyketide|C2",
                "polyketide|C1",
                "polyketide|D2",
                "polyketide|D1",
                "polyketide|C2",
                "polyketide|A1",
                "peptide|pubchem|205",
            ],
            None,
        ),
        (
            "indanomycin",
            r"CC[C@H]1CC[C@@H]2[C@@H]1C=C[C@H]([C@H]2C(=O)C3=CC=CN3)/C=C/C=C(\CC)/[C@H]4[C@H](CC[C@@H](O4)[C@@H](C)C(=O)O)C",
            [
                "polyketide|A1",
                "polyketide|C1",
                "polyketide|D4",
                "polyketide|C1",
                "polyketide|C1",
                "polyketide|C1",
                "polyketide|C4",
                "polyketide|B2",
                "polyketide|D1",
                "polyketide|C2",
            ],
            None,
        ),
        (
            "ircinianin",
            r"C[C@H]1CCC2[C@@H]1C3(C(C=C2C)/C=C(\C)/CCCC4=COC=C4)C(=C(C(=O)O3)C)O",
            [
                "polyketide|D1",
                "polyketide|D2",
                "polyketide|C1",
                "polyketide|C2",
                "polyketide|C1",
                "polyketide|D2",
                "polyketide|D11",
                "polyketide|A2",
            ],
            None,
        ),
        (
            "iriomoteolide_1a",
            r"C[C@H]1C/C=C/[C@@]([C@@]2(CC(=C)C[C@@H](O2)C/C=C/[C@@H]([C@@H](/C(=C\C(=O)O[C@@H]1C[C@H](C)[C@H](C)O)/C)C)O)O)(C)O",
            [
                "polyketide|D8",
                "polyketide|B2",
                "polyketide|C1",
                "polyketide|D7",
                "polyketide|A1",
                "polyketide|A8",
                "polyketide|B1",
                "polyketide|C1",
                "polyketide|B2",
                "polyketide|C6",
            ],
            None,
        ),
        (
            "iriomoteolide_3a",
            r"C/C=C/C/C=C/[C@H](C)C[C@@H]([C@@H]1C[C@H]2[C@@H](O2)/C=C/[C@@H]([C@H](/C=C/C[C@@H](CC(=O)O1)C)O)O)O",
            [
                "polyketide|C1",
                "polyketide|C1",
                "polyketide|D8",
                "polyketide|B5",
                "polyketide|D1",
                "polyketide|C1",
                "polyketide|D11",
                "polyketide|B1",
                "polyketide|C1",
                "polyketide|D8",
            ],
            None,
        ),
        (
            "jerangolid_A",
            r"CC[C@@H]1C(=CC[C@@H](O1)/C(=C/[C@H](C)/C=C/[C@H]2CC(=C(C(=O)O2)CO)OC)/C)C",
            [
                "polyketide|C2",
                "polyketide|C1",
                "polyketide|B2",
                "polyketide|C2",
                "polyketide|C1",
                "polyketide|B1",
                "polyketide|A6",
            ],
            None,
        ),
        (
            "kirromycin",
            r"CC[C@H](C(=O)NC/C=C/C=C(\C)/[C@H]([C@@H](C)[C@H]1[C@H]([C@H]([C@H](O1)/C=C/C=C/C=C(\C)/C(=O)C2=C(C=CNC2=O)O)O)O)OC)[C@@]3([C@@H]([C@@H](C([C@@H](O3)/C=C/C=C/C)(C)C)O)O)O",
            [
                "polyketide|C1",
                "polyketide|C1",
                "polyketide|B3",
                "polyketide|B5",
                "polyketide|A4",
                "peptide|pubchem|750",
                "polyketide|C1",
                "polyketide|C2",
                "polyketide|B2",
                "polyketide|D11",
                "polyketide|B5",
                "polyketide|C1",
                "polyketide|C1",
                "polyketide|C2",
                "polyketide|A1",
                "peptide|pubchem|239",
            ],
            None,
        ),
        (
            "lactimidomycin",
            r"C[C@H]1/C=C\C=C\CC/C=C/C(=O)O[C@H]1/C(=C/[C@H](C)C(=O)C[C@@H](CC2CC(=O)NC(=O)C2)O)/C",
            [
                "polyketide|B1",
                "polyketide|A2",
                "polyketide|C2",
                "polyketide|B2",
                "polyketide|C1",
                "polyketide|C1",
                "polyketide|D1",
                "polyketide|C1",
            ],
            None,
        ),
        (
            "lankamycin",
            r"C[C@@H]1C[C@@H]([C@H]([C@@H](O1)O[C@H]2[C@H](C[C@](C(=O)[C@@H]([C@H]([C@H]([C@H](OC(=O)[C@@H]([C@H]([C@@H]2C)O[C@H]3C[C@@]([C@@H]([C@@H](O3)C)OC(=O)C)(C)OC)C)[C@@H](C)[C@H](C)O)C)OC(=O)C)C)(C)O)C)O)OC",
            [
                "polyketide|B2",
                "polyketide|B2",
                "polyketide|B2",
                "polyketide|A7",
                "polyketide|D2",
                "polyketide|B2",
                "polyketide|B2",
            ],
            None,
        ),
        (
            "latrunculin",
            r"C[C@H]/1CC[C@@H]2C[C@H](C[C@@](O2)([C@@H]3CSC(=O)N3)O)OC(=O)/C=C(\CC/C=C/C=C1)/C",
            [
                "peptide|pubchem|594",
                "polyketide|A1",
                "polyketide|B1",
                "polyketide|B1",
                "polyketide|D2",
                "polyketide|C1",
                "polyketide|C1",
                "polyketide|D1",
                "polyketide|C6",
            ],
            None,
        ),
        # (
        #     "leiodermatolide",
        #     r"CC[C@H]1[C@@H]([C@@](CC(=O)O1)(C/C=C/C=C(\C)/[C@@H]2[C@H](/C=C\C=C/[C@@H]([C@H]([C@H]([C@@H](/C(=C/CCC(=O)O2)/C)C)O)C)OC(=O)N)C)O)C",
        #     [],
        #     None,
        # ),
        # (
        #     "lovastatin",
        #     r"CC[C@H](C)C(=O)O[C@H]1C[C@H](C=C2[C@H]1[C@H]([C@H](C=C2)C)CC[C@@H]3C[C@H](CC(=O)O3)O)C",
        #     [],
        #     None,
        # ),
        (
            "macrolactin_A",
            r"CC1CCC/C=C/C=C/C(CC(C/C=C\C=C\C(C/C=C/C=C\C(=O)O1)O)O)O",
            [
                "polyketide|B1",
                "polyketide|D1",
                "polyketide|C1",
                "polyketide|C1",
                "polyketide|B1",
                "polyketide|B1",
                "polyketide|C1",
                "polyketide|C1",
                "polyketide|B1",
                "polyketide|C1",
                "polyketide|C1",
            ],
            None,
        ),
        (
            "megalomycin_A",
            r"CC[C@@H]1[C@@]([C@@H]([C@H](C(=O)[C@@H](C[C@@]([C@@H]([C@H]([C@@H]([C@H](C(=O)O1)C)O[C@H]2C[C@@]([C@H]([C@@H](O2)C)O)(C)O)C)O[C@H]3[C@@H]([C@H](C[C@H](O3)C)N(C)C)O)(C)O[C@H]4C[C@H]([C@H]([C@@H](O4)C)O)N(C)C)C)C)O)(C)O",
            [
                "polyketide|B7",
                "polyketide|B2",
                "polyketide|A2",
                "polyketide|D7",
                "polyketide|B2",
                "polyketide|B2",
            ],
            None,
        ),
        (
            "micacocidin_A",
            r"CCCCCC1=C(C(=CC=C1)[O-])C2=N[C@H](CS2)[C@@H]3N([C@@H](CS3)[C@@H](C(C)(C)C4=N[C@@](CS4)(C)C(=O)[O-])O)C",
            [
                "polyketide|D1",
                "polyketide|D1",
                "polyketide|C1",
                "polyketide|C1",
                "polyketide|A1",
                "peptide|pubchem|594",
                "peptide|pubchem|594",
                "polyketide|B3",
                "peptide|pubchem|594",
            ],
            None,
        ),
        (
            "migrastatin",
            r"C[C@@H]1/C=C(\[C@H](OC(=O)/C=C/CC/C=C/[C@@H]([C@H]1O)OC)[C@H](C)C(=O)CCCC2CC(=O)NC(=O)C2)/C",
            [
                "polyketide|D1",
                "polyketide|A2",
                "polyketide|B2",
                "polyketide|C2",
                "polyketide|B5",
                "polyketide|C1",
                "polyketide|D1",
                "polyketide|C1",
            ],
            None,
        ),
        (
            "narbonolide",
            r"CC[C@@H]1[C@@H](/C=C/C(=O)[C@@H](C[C@@H]([C@@H]([C@H](C(=O)[C@H](C(=O)O1)C)C)O)C)C)C",
            [
                "polyketide|B2",
                "polyketide|C1",
                "polyketide|A2",
                "polyketide|D2",
                "polyketide|B2",
                "polyketide|A2",
            ],
            None,
        ),
        (
            "pederin",
            r"C[C@H]1[C@H](O[C@](CC1=C)([C@@H](C(=O)N[C@H]([C@@H]2C[C@H](C([C@H](O2)C[C@@H](COC)OC)(C)C)O)OC)O)OC)C",
            [],
            None,
        ),
        (
            "peloruside_A",
            r"CC[C@@H](CO)/C=C(/C)\[C@@H]1C[C@H](C[C@@H](C([C@@]2([C@@H]([C@@H](C[C@@H](O2)C[C@H]([C@@H](C(=O)O1)O)OC)OC)O)O)(C)C)O)OC",
            [
                "polyketide|D6",
                "polyketide|C2",
                "polyketide|B1",
                "polyketide|B1",
                "polyketide|B3",
                "polyketide|A5",
                "polyketide|B1",
                "polyketide|B1",
                "polyketide|B5",
            ],
            None,
        ),
        (
            "penicillin_G",
            r"CC1([C@@H](N2[C@H](S1)[C@@H](C2=O)NC(=O)CC3=CC=CC=C3)C(=O)O)C",
            [],
            None,
        ),
        (
            "periconiasin_I",
            r"C/C/1=C/C[C@@H](CC(=O)[C@]23[C@@H](C1)[C@H](C(=C([C@H]2[C@@H](NC3=O)CC(C)C)C)C)O)O",
            [],
            None,
        ),
        (
            "ratjadon",
            r"C/C=C/[C@H]1[C@H]([C@@H](C[C@H](O1)[C@@H](/C=C/C=C(\C)/C[C@@H](C)/C=C(/C)\C=C\[C@H]2CC=CC(=O)O2)O)O)C",
            [
                "polyketide|C1",
                "polyketide|B2",
                "polyketide|B1",
                "polyketide|D11",
                "polyketide|C1",
                "polyketide|C2",
                "polyketide|D2",
                "polyketide|C2",
                "polyketide|C1",
                "polyketide|B1",
                "polyketide|C1",
            ],
            None,
        ),
        # (
        #     "salinosporamide", 
        #     r"CCCC1C(=O)NC2(C1(OC2=O)C)C(C3CCCC=C3)O", 
        #     [], 
        #     None
        # ),
        (
            "soraphen_A",
            r"C[C@H]1/C=C/[C@H]([C@H](CCCC[C@H](OC(=O)[C@H]([C@@]2([C@@H]([C@H]([C@@H]([C@H]1O2)C)O)OC)O)C)C3=CC=CC=C3)OC)OC",
            [
                "polyketide|B1",
                "polyketide|D1",
                "polyketide|D5",
                "polyketide|B1",
                "polyketide|C2",
                "polyketide|B2",
                "polyketide|B5",
                "polyketide|A2",
            ],
            None,
        ),
        (
            "spiculoic_acid_A",
            r"CC[C@@H]1[C@@H]2[C@@H]([C@H](C1=O)C)C(=C[C@@]([C@@]2(CC)C(=O)O)(CC)/C=C/C3=CC=CC=C3)CC",
            [
                "polyketide|C4",
                "polyketide|C4",
                "polyketide|C2",
                "polyketide|A4",
                "polyketide|C4",
            ],
            None,
        ),
        (
            "spongidepsin",
            r"C[C@H]1CCC(CC(OC(=O)[C@@H](N(C(=O)[C@H](C1)C)C)CC2=CC=CC=C2)CCCC#C)C",
            [
                "peptide|pubchem|*",
                "polyketide|D1",
                "polyketide|B1",
                "polyketide|D8",
                "polyketide|D2",
                "polyketide|D2",
                "peptide|pubchem|994",
            ],
            None,
        ),
        (
            "thailanstatin_A",
            r"C[C@H]1C[C@H]([C@H](O[C@H]1C/C=C(\C)/C=C/[C@@H]2[C@H]([C@@]3(C[C@H](O2)CC(=O)O)CO3)O)C)NC(=O)/C=C\[C@H](C)OC(=O)C",
            [
                "polyketide|C1",
                "peptide|pubchem|205",
                "polyketide|D2",
                "polyketide|C1",
                "polyketide|C2",
                "polyketide|C1",
                "polyketide|B5",
                "polyketide|A8",
                "polyketide|C1",
            ],
            None,
        ),
        # (
        #     "theopederin_A",
        #     r"C[C@H]1[C@H](O[C@](CC1=C)([C@H](C(=O)N[C@@H]2[C@@H]3[C@@H]([C@H](C([C@H](O3)C[C@H]4CCCC(O4)O)(C)C)OC)OCO2)O)OC)C",
        #     [],
        #     None,
        # ),
        # (
        #     "theopederin_B",
        #     r"C[C@H]1[C@H](O[C@](CC1=C)([C@H](C(=O)N[C@@H]2[C@@H]3[C@@H]([C@H](C([C@H](O3)C[C@@H](CCCC(=O)OC)O)(C)C)OC)OCO2)O)OC)C",
        #     [],
        #     None,
        # ),
        (
            "thermolide_A",
            r"C[C@@H]1C[C@H]([C@@H](OC(=O)[C@H](NC(=O)C[C@H](C[C@@H]1O)O)C)[C@@H](C)C[C@H](C)[C@@H]([C@H](C)[C@@H](C[C@H](C)O)OC(=O)C)O)C",
            [
                "polyketide|B1",
                "polyketide|B2",
                "polyketide|B2",
                "polyketide|D2",
                "polyketide|B2",
                "polyketide|D2",
                "polyketide|B1",
                "polyketide|B1",
                "peptide|pubchem|602",
            ],
            None,
        ),
        # (
        #     "thiocoraline",
        #     r"CN1C2CSSCC(C(=O)N(C(C(=O)SCC(C(=O)NCC1=O)NC(=O)C3=NC4=CC=CC=C4C=C3O)CSC)C)N(C(=O)CNC(=O)C(CSC(=O)C(N(C2=O)C)CSC)NC(=O)C5=NC6=CC=CC=C6C=C5O)C",
        #     [],
        #     RetroMolTimeoutError
        # ),
        (
            "zincophorin",
            r"CCC[C@@H](C)/C=C(\C)/[C@@H]([C@H](C)/C=C/CC[C@H]([C@H](C)[C@@H]([C@H](C)[C@@H]([C@H](C)[C@@H]1[C@H](CC[C@H](O1)[C@H](C)C(=O)O)C)O)O)O)O",
            [
                "polyketide|D2",
                "polyketide|C2",
                "polyketide|B2",
                "polyketide|C1",
                "polyketide|D1",
                "polyketide|B2",
                "polyketide|B2",
                "polyketide|B2",
                "polyketide|C2",
                "polyketide|D1",
                "polyketide|B2",
            ],
            None,
        ),
        (
            "zwittermicin_A",
            r"C([C@H]([C@H](CO)N)O)[C@H]([C@H]([C@H]([C@@H](C(=O)N[C@@H](CNC(=O)N)C(=O)N)O)O)N)O",
            [
                "peptide|pubchem|617",
                "polyketide|B1",
                "polyketide|B12",
                "polyketide|B5",
                "peptide|pubchem|364",
            ],
            None,
        ),
        # (
        #     "enterobactin",
        #     r"C1C(C(=O)OCC(C(=O)OCC(C(=O)O1)NC(=O)C2=C(C(=CC=C2)O)O)NC(=O)C3=C(C(=CC=C3)O)O)NC(=O)C4=C(C(=CC=C4)O)O",
        #     [
        #         [],
        #         [],
        #         [],
        #     ],
        #     None,
        # )
    ],
)
def test_modular_natural_produts(
    compound_name: str, test: str, expected: ty.List[str], expected_raises: ty.Any
):
    """Integration test for parsing and sequencing a natural product compound.

    :param test: The test input.
    :type test: str
    :param expected: The expected result.
    :type expected: List[str]
    :param expected_raises: The expected exception.
    :type expected_raises: Any
    """
    molecule = Molecule("test_input", test)

    result = parse_mol_timed(molecule)

    if expected_raises is not None:
        with pytest.raises(expected_raises):
            assert not result.success
    else:
        assert result.success

        sequences = result.sequences

        if len(sequences) > 0:
            motif_codes = [seq["motif_code"] for seq in sequences]
            print(motif_codes)
            assert expected in motif_codes

        else:
            assert sequences == []
