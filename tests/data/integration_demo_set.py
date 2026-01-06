# -*- coding: utf-8 -*-

"""Integration test cases for molecule coverage scoring and optimal mapping identification."""

# Cases are formatted as:
#
# (
#    name,
#    smiles,
#    expected coverage_score,
#    list of found monomers
# ),
#
CASES: list[tuple[str, str, float, list[str]]] = [
    (
        "10-deoxymethynolide",
        r"CC[C@@H]1[C@@H](/C=C/C(=O)[C@@H](C[C@@H]([C@@H]([C@H](C(=O)O1)C)O)C)C)C",
        1.0,
        ["A2", "B2", "B2", "C1", "D2", "propanoic acid"],
    ),
    (
        "13-deoxytedanolide",
        r"C/C=C\[C@H](C)[C@@H]1[C@](O1)(C)[C@H]([C@H]2COC(=O)[C@@H]([C@H]([C@@H](C(=O)[C@@H]([C@H](/C(=C/[C@@H](C(=O)CC[C@H](C2=O)C)C)/C)O)C)C)OC)O)O",
        1.0,
        ["A2", "A2", "A2", "B2", "B5", "B7", "C2", "C2", "C2", "D1", "methylation", "oxidation", "propanoic acid"],
    ),
    (
        "2-deoxystreptamine",
        r"C1C(C(C(C(C1N)O)O)O)N",
        1.0,
        ["4,6-diaminocyclohexane-1,2,3-triol"],
    ),
    (
        "6-deoxyerythronolide B",
        r"CC[C@@H]1[C@@H]([C@@H]([C@H](C(=O)[C@@H](C[C@@H]([C@@H]([C@H]([C@@H]([C@H](C(=O)O1)C)O)C)O)C)C)C)O)C",
        1.0,
        ["A2", "B2", "B2", "B2", "B2", "D2", "propanoic acid"],
    ),
    (
        "AF-toxin",
        r"CCC(C)C(C(=O)OC(/C=C/C=C/C=C/C(=O)O)C1(CO1)C)OC(=O)C(C(C)(C)O)O",
        1.0,
        ["acetic acid", "B10", "acetic acid", "D17", "2-hydroxy-3-methylpentanoic acid", "C1", "C1", "C1", "oxidation"]
    ),
    (
        "abyssomicin C",
        r"C[C@@H]1C[C@]23OC(=O)C4=C2OC1[C@H](O)C3\C=C\C(=O)[C@@H](C)C[C@@H](C)C4=O",
        0.72,
        ["A2", "C1", "C1", "C1", "D2", "acetic acid", "glyceric acid", "oxidation"],
    ),
    (
        "atrop-abyssocymin C",
        r"CC1CC23OC(=O)C4=C2OC1C(O)C3\C=C/C(=O)C(C)CC(C)C4=O",
        0.72,
        ["A2", "C1", "C1", "C1", "D2", "acetic acid", "glyceric acid", "oxidation"],
    ),
    (
        "aculeximycin",
        r"CCCC(O[C@H]1C[C@](C)(N)[C@H](O)[C@H](C)O1)C(C)C(O)C(CC)\C=C\C(O)C(C)C1C\C=C(C)\C(O)C(C)C(CC(O)C(C)C(O)CC2CC(O)C(O)C(O)(CC(O[C@@H]3O[C@H](C)[C@@H](O)[C@H](O[C@H]4C[C@@H](N)[C@H](O)[C@@H](C)O4)[C@H]3O[C@@H]3O[C@H](C)[C@@H](O)[C@H](O)[C@H]3O)C(C)CCC(O)CC(O)C\C=C(CC)\C(=O)O1)O2)O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]1O",
        1.0,
        ["4-amino-4,6-dimethyloxane-2,5-diol", "A1", "B1", "B1", "B1", "B1", "B1", "B1", "B2", "B2", "B2", "B2", "B2", "B4", "B5", "C1", "C2", "C4", "D1", "D1", "acetic acid", "butanoic acid", "glucose", "rhamnose", "rhamnose", "sugar"],
    ),
    (
        "acutiphycin",
        r"CCCCC[C@@H]1C/C=C(\[C@H](C(C(=O)[C@H](/C=C(\[C@@H]2C[C@@H](C[C@@](O2)(CC(=O)O1)O)O)/C)C)(C)C)O)/C",
        1.0,
        ["A1", "A2", "B1", "B1", "B1", "B3", "C2", "C2", "D1", "D1", "acetic acid", "butanoic acid", "hexanoic acid"],
    ),
    (
        "aflatoxin G1",
        r"COC1=C2C3=C(C(=O)OCC3)C(=O)OC2=C4[C@H]5C=CO[C@H]5OC4=C1",
        0.04,
        ["methylation"],
    ),
    (
        "alternapyrone",
        r"CCC(C)/C=C(\C)/CCCC(C)/C=C(\C)/C=C(\C)/CC(C)C1=C(C(=C(C(=O)O1)C)O)C",
        1.0,
        ["A2", "A2", "C2", "C2", "C2", "D1", "D2", "D2", "D2", "acetic acid"],
    ),
    (
        "amicoumacin",
        r"CC(C)C[C@@H]([C@@H]1CC2=C(C(=CC=C2)O)C(=O)O1)NC(=O)[C@H]([C@H]([C@H](CC(=O)N)N)O)O",
        1.0,
        ["A1", "B1", "B5", "C1", "C1", "asparagine", "leucine"],
    ),
    (
        "amphidinolide J",
        r"CCC/C=C/[C@@H](C)[C@H]1C(/C=C\C([C@H](C=CCCC(=C)[C@H](CC(=O)O1)C)O)C)O",
        1.0,
        ["B1", "B1", "C1", "C1", "C2", "D1", "D10", "D15", "D8", "acetic acid", "trans-2-hexanoic acid", "butanoic acid"],
    ),
    (
        "amphidinolide P",
        r"C[C@@H]1C(=C)C[C@H]2[C@H]3[C@@H](O3)CC(=C)/C=C/[C@H](OC(=O)C[C@@]1(O2)O)[C@H](C)C(=C)C",
        1.0,
        ["2-methylprop-2-enoic acid", "A1", "A8", "A9", "B1", "C1", "C1", "D10", "oxidation"],
    ),
    (
        "ansamitocin P-3 ",
        r"C[C@@H]1[C@@H]2C[C@]([C@@H](/C=C/C=C(/CC3=CC(=C(C(=C3)OC)Cl)N(C(=O)C[C@@H]([C@]4([C@H]1O4)C)OC(=O)C(C)C)C)\C)OC)(NC(=O)O2)O",
        1.0,
        ["3-amino-5-hydroxybenzoic acid", "A1", "B1", "B2", "C1", "C2", "D11", "D2", "carbamic acid", "chlorination", "isobutyric acid", "methylation", "methylation", "methylation", "oxidation"],
    ),
    (
        "anthracimycin",
        r"C[C@@H]1/C=C\C=C\[C@H](OC(=O)[C@@H](C(=O)/C=C(/[C@H]2[C@@H]1C=C[C@@H]3[C@@H]2CC=C(C3)C)\O)C)C",
        1.0,
        ["A1", "A2", "B1", "C1", "C1", "C1", "C1", "C1", "C2", "D2", "acetic acid"],
    ),
    (
        "apoptolidin",
        r"COC[C@@H](C[C@H]1O[C@@](O)([C@H](O)[C@@H]2C[C@H](OC)[C@@H](O)CC\C=C(/C)\C=C\[C@@H](O[C@@H]3O[C@@H](C)[C@H](OC)[C@@H](O)[C@@H]3O)[C@H](C)\C=C(/C)\C=C(/C)\C=C(C)\C(=O)O2)[C@H](C)[C@@H](O)[C@H]1C)O[C@H]1C[C@](C)(O)[C@@H](O[C@H]2C[C@@H](OC)[C@H](O)[C@@H](C)O2)[C@H](C)O1",
        1.0,
        ["4,6-dimethyloxane-2,4,5-triol", "6-methyloxane-2,4,5-triol", "A2", "B1", "B1", "B2", "B2", "B5", "C1", "C2", "C2", "C2", "C2", "D5", "malonic acid", "methanol", "methylation", "methylation", "methylation", "methylation", "rhamnose"],
    ),
    (
        "avilamycin A",
        r"C[C@@H]1[C@H]([C@@H](C[C@@H](O1)O[C@@H]2[C@H](OC3(C[C@@H]2O)O[C@@H]4[C@H](O[C@H](C[C@]4(O3)C)O[C@@H]5[C@H]([C@@H](O[C@@H]([C@@H]5OC)C)O[C@@H]6[C@H](O[C@H]([C@H]([C@H]6O)OC)O[C@H]7[C@@H]([C@H]8[C@H](CO7)O[C@@]9(O8)C1C([C@@]([C@H](O9)C)(C(=O)C)O)OCO1)OC(=O)C(C)C)COC)O)C)C)O)OC(=O)C1=C(C(=C(C(=C1OC)Cl)O)Cl)C",
        1.0,
        ["4,6-dimethyloxane-2,4,5-triol", "6-methyloxane-2,4,5-triol", "6-methyloxane-2,4,5-triol", "A1", "A1", "C1", "acetic acid", "arabinose", "chlorination", "chlorination", "glucose", "isobutyric acid", "methanol", "methylation", "methylation", "methylation", "methylation", "orsellinic acid", "rhamnose", "sugar"],
    ),
    (
        "avilamycin C",
        r"C[C@@H]1[C@H]([C@@H](C[C@@H](O1)O[C@@H]2[C@H](OC3(C[C@H]2O)O[C@@H]4[C@H](O[C@H](C[C@]4(O3)C)O[C@@H]5[C@H]([C@@H](O[C@@H]([C@@H]5OC)C)O[C@@H]6[C@H](O[C@H]([C@H]([C@H]6O)OC)O[C@H]7[C@@H]([C@H]8[C@H](CO7)O[C@@]9(O8)[C@H]1[C@H]([C@@]([C@H](O9)C)(C(C)O)O)OCO1)OC(=O)C(C)C)COC)O)C)C)O)OC(=O)C1=C(C(=C(C(=C1OC)Cl)O)Cl)C",
        1.0,
        ["4,6-dimethyloxane-2,4,5-triol", "6-methyloxane-2,4,5-triol", "6-methyloxane-2,4,5-triol", "A1", "A1", "C1", "acetic acid", "arabinose", "chlorination", "chlorination", "glucose", "isobutyric acid", "methanol", "methylation", "methylation", "methylation", "methylation", "orsellinic acid", "rhamnose", "sugar"],
    ),
    (
        "bitungolide F",
        r"CC[C@@H]1C=CC(=O)O[C@@H]1[C@H](C)CC[C@H](C[C@@H](/C=C/C=C/C2=CC=CC=C2)O)O",
        1.0,
        ["B1", "B1", "B4", "C1", "C1", "D2", "cinnamic acid"],
    ),
    (
        "borrelidin",
        r"C[C@H]1C[C@H](C[C@@H]([C@H](/C(=C\C=C\C[C@H](OC(=O)C[C@@H]([C@H](C1)C)O)[C@@H]2CCC[C@H]2C(=O)O)/C#N)O)C)C",
        1.0,
        ["B1", "B1", "B2", "C1", "C1", "D2", "D2", "D2", "cyanide", "cyclopentane-1,2-dicarboxylic acid"],
    ),
    (
        "butirosin A",
        r"C1[C@@H]([C@H]([C@@H]([C@H]([C@@H]1NC(=O)[C@H](CCN)O)O)O[C@H]2[C@@H]([C@H]([C@H](O2)CO)O)O)O[C@@H]3[C@@H]([C@H]([C@@H]([C@H](O3)CN)O)O)N)N",
        1.0,
        ["3-amino-6-(aminomethyl)oxane-2,4,5-triol", "4,6-diaminocyclohexane-1,2,3-triol", "4-amino-2-hydroxybutanoic acid", "D5", "glycine", "ribose"],
    ),
    (
        "butirosin B",
        r"C1[C@@H]([C@H]([C@@H]([C@H]([C@@H]1NC(=O)[C@H](CCN)O)O)O[C@H]2[C@@H]([C@@H]([C@H](O2)CO)O)O)O[C@@H]3[C@@H]([C@H]([C@@H]([C@H](O3)CN)O)O)N)N",
        1.0,
        ["3-amino-6-(aminomethyl)oxane-2,4,5-triol", "4,6-diaminocyclohexane-1,2,3-triol", "4-amino-2-hydroxybutanoic acid", "D5", "glycine", "ribose"],
    ),
    (
        "calicheamicin",
        r"CCN[C@H]1CO[C@H](C[C@@H]1OC)O[C@@H]2[C@H]([C@@H]([C@H](OC2O[C@H]3C#C/C=C\C#C[C@]\4(CC(=O)C(=C3/C4=C\CSSSC)NC(=O)OC)O)C)NO[C@H]5C[C@@H]([C@@H]([C@H](O5)C)SC(=O)C6=C(C(=C(C(=C6OC)OC)O[C@H]7[C@@H]([C@@H]([C@H]([C@@H](O7)C)O)OC)O)I)C)O)O",
        0.74,
        # Default rule set is not able to parse the enediyne core, but should at least identify the sugar parts
        ["5-amino-6-methyloxane-2,3,4-triol", "6-methyl-5-sulfanyloxane-2,4-diol", "A1", "A5", "C1", "acetic acid", "carbonic acid", "ethanol", "iodination", "methylation", "methylation", "methylation", "methylation", "methylation", "methylation", "rhamnose", "sugar"]
    ),
    (
        "callystatin",
        r"CC[C@H](C)[C@H]([C@H](C)C(=O)[C@H](C)/C=C(\C)/C=C/C[C@@H](C)/C=C(/CC)\C=C\[C@H]1CC=CC(=O)O1)O",
        1.0,
        ["A2", "B1", "B2", "C1", "C1", "C1", "C2", "C4", "D2", "D2", "acetic acid"],
    ),
    (
        "carolacton",
        r"C[C@@H]\1CCC[C@@H]([C@H](OC(=O)[C@@H]([C@@H](/C=C1)O)O)/C(=C/[C@@H](C)C(=O)[C@H](C)[C@@H](CC(=O)O)OC)/C)C",
        1.0,
        # This parsing is ambiguous... can parse from two sides and neither is correct/wrong
        ["A2", "B1", "B2", "B5", "C1", "C2", "D2", "D2", "malonic acid", "methylation"],
    ),
    (
        "chaetoglobosin A",
        r"C[C@H]\1C/C=C/[C@H]2[C@H]3[C@](O3)([C@H]([C@@H]4[C@@]2(C(=O)/C=C/C(=O)[C@@H](/C(=C1)/C)O)C(=O)N[C@H]4CC5=CNC6=CC=CC=C65)C)C",
        1.0,
        ["A1", "B11", "C1", "C1", "C1", "C2", "C2", "D2", "acetic acid", "oxidation", "tryptophan"],
    ),
    (
        "chichorine",
        r"CC1=C(O)C=C2C(CNC2=O)=C1OC",
        1.0,
        ["C1", "D11", "D17", "glycine", "methylation"],
    ),
    (
        "chlorotonil A",
        r"C[C@@H]1/C=C\C=C\[C@@H](OC(=O)[C@H](C(=O)C(C(=O)[C@@H]2[C@H]1C=C[C@H]3[C@H]2[C@@H](C=C(C3)C)C)(Cl)Cl)C)C",
        1.0,
        ["A1", "A2", "B1", "C1", "C1", "C1", "C1", "C2", "C2", "D2", "acetic acid", "chlorination", "chlorination"],
    ),
    (
        "chlorothricin",
        r"C[C@@H](C(C(O)=O)=C1)C[C@@]2(C(O)=C(O3)C(O2)=O)[C@@H]1C=CCCCC[C@@H]4C=C[C@@]([C@@H](O[C@@H]5C[C@H](O)[C@@H](O[C@H]6O[C@@H](C)[C@H](O)[C@@H](OC(C7=C(C)C(Cl)=CC=C7OC)=O)C6)[C@@H](C)O5)CCC8)([H])[C@]8([H])[C@@H]4C3=O",
        1.0,
        ["6-methyloxane-2,4,5-triol", "6-methyloxane-2,4,5-triol", "6-methylsalicylic acid", "A1", "B1", "C1", "C1", "C1", "C1", "C1", "C1", "C1", "C13", "D1", "D1", "D1", "acetic acid", "acetic acid", "chlorination", "glyceric acid", "glycolic acid", "methylation"],
    ),
    (
        "coelimycin P1",
        r"C/C=C/C(=O)/C/1=C/C(=C\2/C=CCCN2)/SC[CH](C(=O)O1)NC(=O)C",
        0.38,
        ["acetic acid", "cysteine"],
    ),
    (
        "compactin",
        r"CC[C@H](C)C(=O)O[C@H]1CCC=C2C=C[C@H](C)[C@H](CC[C@@H](O)C[C@@H](O)CC(O)=O)[C@@H]12",
        1.0,
        ["B1", "B1", "C1", "C1", "C1", "C1", "D1", "D2", "D5", "acetic acid", "acetic acid"],
    ),
    (
        "cremimycin",
        r"CCCCCCC1CC(CCCC(C/C(O)=C2C(/C=C(/C=C\C=C/C(N1)=O)C)CC(C\2=O)O[C@H]3C[C@@H]([C@@H]([C@H](O3)C)O)OC)O)=O",
        1.0,
        ["6-methyloxane-2,4,5-triol", "A", "A1", "A1", "B1", "B1", "C1", "C1", "C1", "C2", "D1", "D1", "D1", "heptanoic acid", "methylation", "oxidation", "pentanoic acid", "propanoic acid"],
    ),
    (
        "deschlorothricin",
        r"C[C@@H](C(C(O)=O)=C1)C[C@@]2(C(O)=C(O3)C(O2)=O)[C@@H]1C=CCCCC[C@@H]4C=C[C@@]([C@@H](O[C@@H]5C[C@H](O)[C@@H](O[C@H]6O[C@@H](C)[C@H](O)[C@@H](OC(C7=C(C)C=CC=C7OC)=O)C6)[C@@H](C)O5)CCC8)([H])[C@]8([H])[C@@H]4C3=O",
        1.0,
        ["6-methyloxane-2,4,5-triol", "6-methyloxane-2,4,5-triol", "6-methylsalicylic acid", "A1", "B1", "C1", "C1", "C1", "C1", "C1", "C1", "C1", "C13", "D1", "D1", "D1", "acetic acid", "acetic acid", "glyceric acid", "glycolic acid", "methylation"],
    ),
    (
        "daptomycin",
        r"CCCCCCCCCC(=O)N[C@@H](CC1=CNC2=CC=CC=C21)C(=O)N[C@H](CC(=O)N)C(=O)N[C@@H](CC(=O)O)C(=O)N[C@H]3[C@H](OC(=O)[C@@H](NC(=O)[C@@H](NC(=O)[C@H](NC(=O)CNC(=O)[C@@H](NC(=O)[C@H](NC(=O)[C@@H](NC(=O)[C@@H](NC(=O)CNC3=O)CCCN)CC(=O)O)C)CC(=O)O)CO)[C@H](C)CC(=O)O)CC(=O)C4=CC=CC=C4N)C",
        1.0,
        ["3-methylglutamic acid", "D1", "D1", "D1", "D1", "acetic acid", "alanine", "asparagine", "aspartic acid", "aspartic acid", "aspartic acid", "butanoic acid", "decanoic acid", "glycine", "glycine", "hexanoic acid", "kynurenine", "octanoic acid", "ornithine", "serine", "threonine", "tryptophan"]
    ),
    (
        "dictyostatin",
        r"C[C@H]1CC[C@H]([C@@H]([C@@H](OC(=O)/C=C\C=C\[C@H]([C@H](C[C@@H](/C=C\[C@@H]([C@@H]([C@H](C1)C)O)C)O)O)C)[C@@H](C)/C=C\C=C)C)O",
        1.0,
        ["B1", "B1", "B2", "B2", "B2", "C1", "C1", "C1", "C1", "C2", "D2", "D2", "acetic acid"],
    ),
    (
        "discodermolide",
        r"C[C@H]1[C@@H](OC(=O)[C@@H]([C@H]1O)C)C[C@@H](/C=C\[C@H](C)[C@@H]([C@@H](C)/C=C(/C)\C[C@H](C)[C@H]([C@H](C)[C@H]([C@@H](C)/C=C\C=C)OC(=O)N)O)O)O",
        1.0,
        ["B1", "B2", "B2", "B2", "B2", "B2", "C1", "C1", "C2", "C2", "D2", "acetic acid", "carbamic acid"],
    ),
    (
        "epothilone",
        r"C[C@H]1CCC[C@@H]2[C@@H](O2)C[C@H](OC(=O)C[C@H](C(C(=O)[C@@H]([C@H]1O)C)(C)C)O)/C(=C/C3=CSC(=N3)C)/C",
        1.0,
        ["A3", "B1", "B1", "B2", "C1", "C2", "D1", "D2", "acetic acid", "cysteine", "oxidation"],
    ),
    (
        "erythromycin",
        r"CC[C@@H]1[C@@]([C@@H]([C@H](C(=O)[C@@H](C[C@@]([C@@H]([C@H]([C@@H]([C@H](C(=O)O1)C)O[C@H]2C[C@@]([C@H]([C@@H](O2)C)O)(C)OC)C)O[C@H]3[C@@H]([C@H](C[C@H](O3)C)N(C)C)O)(C)O)C)C)O)(C)O",
        1.0,
        ["4,6-dimethyloxane-2,4,5-triol", "4-dimethylamino-6-methyloxane-2,3-diol", "A2", "B2", "B2", "B2", "B6", "D6", "methylation", "propanoic acid"]
    ),
    (
        "georatusin",
        r"CC[C@@H](C)C[C@H](C)[C@@H]1[C@H](C[C@@H]([C@@H]2[C@H](C[C@@H]([C@@](O2)([C@H](C(=O)N[C@@H](C(=O)O1)CC3=CNC4=CC=CC=C43)C)O)C)C)C)C",
        1.0,
        ["A2", "B2", "B2", "D2", "D2", "D2", "D2", "acetic acid", "tryptophan"],
    ),
    (
        "gephyronic acid",
        r"C[C@@H]1[C@@H](O[C@@](C([C@H]1OC)(C)C)([C@@H](C)C[C@H](C)[C@@H]([C@@]2([C@H](O2)[C@@H](C)C=C(C)C)C)O)O)CC(=O)O",
        1.0,
        ["A3", "B1", "B2", "B2", "C2", "C2", "D2", "isobutyric acid", "methylation", "oxidation"],
    ),
    (
        "harzianic acid",
        r"CCC/C=C/C=C/C(=C\1/C(=O)C(N(C1=O)C)CC(C(C)C)(C(=O)O)O)/O",
        1.0,
        ["A1", "C1", "C1", "D1", "acetic acid", "artificial amino acid harzianic acid", "butanoic acid", "methylation", "trans-2-hexanoic acid"]
    ),
    (
        "herboxidiene",
        r"C[C@H]1CC[C@@H](O[C@@H]1/C(=C/C=C/[C@@H](C)C[C@@]2([C@H](O2)[C@H](C)[C@H]([C@@H](C)O)OC)C)/C)CC(=O)O",
        1.0,
        ["B2", "B2", "C1", "C1", "C2", "C2", "D1", "D2", "lactic acid", "methylation", "oxidation"],
    ),
    (
        "hydroxystreptomycin",
        r"CN[C@H]1[C@H](O)[C@@H](O)[C@H](CO)O[C@H]1O[C@H]1[C@H](O[C@H]2[C@H](O)[C@@H](O)[C@H](NC(N)=N)[C@@H](O)[C@@H]2NC(N)=N)O[C@@H](CO)[C@]1(O)C=O",
        1.0,
        ["2-[3-(diaminomethylideneamino)-2,4,5,6-tetrahydroxycyclohexyl]guanidine", "glucosamine", "methylation", "sugar"],
    ),
    (
        "hymenosetin",
        r"C/C=C/[C@@H]1C(=C[C@@H]2C[C@@H](CC[C@H]2[C@]1(C)/C(=C\3/C(=O)[C@@H](NC3=O)[C@@H](C)O)/O)C)C",
        1.0,
        ["A1", "C1", "C1", "C2", "C2", "D1", "D2", "acetic acid", "threonine"],
    ),
    (
        "indanomycin",
        r"CC[C@H]1CC[C@@H]2[C@@H]1C=C[C@H]([C@H]2C(=O)C3=CC=CN3)/C=C/C=C(\CC)/[C@H]4[C@H](CC[C@@H](O4)[C@@H](C)C(=O)O)C",
        1.0,
        ["A1", "B2", "C1", "C1", "C1", "C1", "C2", "C4", "D1", "D4", "pyrrole-2-carboxylic acid"],
    ),
    (
        "ircinianin",
        r"C[C@H]1CCC2[C@@H]1C3(C(C=C2C)/C=C(\C)/CCCC4=COC=C4)C(=C(C(=O)O3)C)O",
        1.0,
        ["A2", "C1", "C1", "C2", "D1", "D11", "D2", "D2", "furan-3-carboxylic acid"],
    ),
    (
        "iriomoteolide 1a",
        r"C[C@H]1C/C=C/[C@@]([C@@]2(CC(=C)C[C@@H](O2)C/C=C/[C@@H]([C@@H](/C(=C\C(=O)O[C@@H]1C[C@H](C)[C@H](C)O)/C)C)O)O)(C)O",
        1.0,
        ["A1", "A8", "A8", "B1", "B2", "B2", "C1", "C1", "D6", "D8", "lactic acid"],
    ),
    (
        "iriomoteolide 3a",
        r"C[C@H]1C/C=C/[C@@]([C@@]2(CC(=C)C[C@@H](O2)C/C=C/[C@@H]([C@@H](/C(=C\C(=O)O[C@@H]1C[C@H](C)[C@H](C)O)/C)C)O)O)(C)O",
        1.0,
        ["A1", "A8", "A8", "B1", "B2", "B2", "C1", "C1", "D6", "D8", "lactic acid"],
    ),
    (
        "jerangolid A",
        r"CC[C@@H]1C(=CC[C@@H](O1)/C(=C/[C@H](C)/C=C/[C@H]2CC(=C(C(=O)O2)CO)OC)/C)C",
        1.0,
        ["A7", "B1", "B2", "C1", "C2", "C2", "D1", "methylation", "propanoic acid"],
    ),
    (
        "kirromycin",
        r"CC[C@H](C(=O)NC/C=C/C=C(\C)/[C@H]([C@@H](C)[C@H]1[C@H]([C@H]([C@H](O1)/C=C/C=C/C=C(\C)/C(=O)C2=C(C=CNC2=O)O)O)O)OC)[C@@]3([C@@H]([C@@H](C([C@@H](O3)/C=C/C=C/C)(C)C)O)O)O",
        1.0,
        ["A1","A4","B2","B3","B5","B5","C1","C1","C1","C1","C1","C2","C2","D11","acetic acid","beta-alanine","glycine","methylation"],
    ),
    (
        "lactimidomycin",
        r"C[C@H]1/C=C\C=C\CC/C=C/C(=O)O[C@H]1/C(=C/[C@H](C)C(=O)C[C@@H](CC2CC(=O)NC(=O)C2)O)/C",
        1.0,
        ["2-(2,6-dioxopiperidin-4-yl)acetic acid", "A2", "B1", "B2", "C1", "C1", "C1", "C2", "D1"],
    ),
    (
        "lankamycin",
        r"C[C@@H]1C[C@@H]([C@H]([C@@H](O1)O[C@H]2[C@H](C[C@](C(=O)[C@@H]([C@H]([C@H]([C@H](OC(=O)[C@@H]([C@H]([C@@H]2C)O[C@H]3C[C@@]([C@@H]([C@@H](O3)C)OC(=O)C)(C)OC)C)[C@@H](C)[C@H](C)O)C)OC(=O)C)C)(C)O)C)O)OC",
        1.0,
        ["4,6-dimethyloxane-2,4,5-triol", "6-methyloxane-2,3,4-triol", "A6", "B2", "B2", "B2", "B2", "B2", "D2", "acetic acid", "acetic acid", "acetic acid", "methylation", "methylation"],
    ),
    (
        "latrunculin",
        r"C[C@H]1CC[C@@H]2C[C@H](C[C@@](O2)([C@@H]3CSC(=O)N3)O)OC(=O)C=C(CCC=CC=C1)C",
        1.0,
        ["A1", "A8", "B1", "B1", "C1", "C1", "D1", "D2", "carbonic acid", "cysteine"],
    ),
    (
        "leiodermatolide",
        r"CC[C@H]1[C@@H]([C@@](CC(=O)O1)(C/C=C/C=C(\C)/[C@@H]2[C@H](/C=C\C=C/[C@@H]([C@H]([C@H]([C@@H](/C(=C/CCC(=O)O2)/C)C)O)C)OC(=O)N)C)O)C",
        0.88,
        ["A8", "B2", "B2", "B2", "B2", "C1", "C1", "C1", "C2", "D1", "carbamic acid", "malonic acid", "propanoic acid"],
    ),
    (
        "lovastatin",
        r"CC[C@H](C)C(=O)O[C@H]1C[C@H](C=C2[C@H]1[C@H]([C@H](C=C2)C)CC[C@@H]3C[C@H](CC(=O)O3)O)C",
        1.0,
        ["B1", "B1", "C1", "C1", "C1", "C2", "D1", "D2", "D5", "acetic acid", "acetic acid"],
    ),
    (
        "macrolactin A",
        r"CC1CCCC=CC=CC(CC(CC=CC=CC(CC=CC=CC(=O)O1)O)O)O",
        1.0,
        ["B1", "B1", "B1", "B1", "C1", "C1", "C1", "C1", "C1", "C1", "D1", "acetic acid"],
    ),
    (
        "maytansine",
        r"C[C@@H]1[C@@H]2C[C@]([C@@H](/C=C/C=C(/CC3=CC(=C(C(=C3)OC)Cl)N(C(=O)C[C@@H]([C@]4([C@H]1O4)C)OC(=O)[C@H](C)N(C)C(=O)C)C)\C)OC)(NC(=O)O2)O",
        1.0,
        ["3-amino-5-hydroxybenzoic acid", "A1", "B1", "B2", "C1", "C2", "D11", "D2", "acetic acid", "alanine", "carbamic acid", "chlorination", "methylation", "methylation", "methylation", "methylation", "oxidation"],
    ),
    (
        "megalomycin A",
        r"CC[C@@H]1[C@@]([C@@H]([C@H](C(=O)[C@@H](C[C@@]([C@@H]([C@H]([C@@H]([C@H](C(=O)O1)C)O[C@H]2C[C@@]([C@H]([C@@H](O2)C)O)(C)O)C)O[C@H]3[C@@H]([C@H](C[C@H](O3)C)N(C)C)O)(C)O[C@H]4C[C@H]([C@H]([C@@H](O4)C)O)N(C)C)C)C)O)(C)O",
        1.0,
        ["4,6-dimethyloxane-2,4,5-triol", "4-dimethylamino-6-methyloxane-2,3-diol", "A2", "B2", "B2", "B2", "B6", "D6", "methylation", "methylation", "propanoic acid", "sugar"],
    ),
    (
        "micacocidin A",
        r"CCCCCC1=C(C(=CC=C1)[O-])C2=N[C@H](CS2)[C@@H]3N([C@@H](CS3)[C@@H](C(C)(C)C4=N[C@@](CS4)(C)C(=O)[O-])O)C",
        1.0,
        ["A1", "B3", "C1", "C1", "D1", "D1", "acetic acid", "cysteine", "cysteine", "cysteine", "methylation", "methylation", "butanoic acid", "hexanoic acid"],
    ),
    (
        "migrastatin",
        r"C[C@@H]1/C=C(\[C@H](OC(=O)/C=C/CC/C=C/[C@@H]([C@H]1O)OC)[C@H](C)C(=O)CCCC2CC(=O)NC(=O)C2)/C",
        1.00,
        ["2-(2,6-dioxopiperidin-4-yl)acetic acid", "A2", "B2", "B5", "C1", "C1", "C2", "D1", "D1", "methylation"],
    ),
    (
        "narbonolide",
        r"CC[C@@H]1[C@@H](/C=C/C(=O)[C@@H](C[C@@H]([C@@H]([C@H](C(=O)[C@H](C(=O)O1)C)C)O)C)C)C",
        1.0,
        ["A2", "A2", "B2", "B2", "C1", "D2", "propanoic acid"],
    ),
    (
        "pederin",
        r"C[C@H]1[C@H](O[C@](CC1=C)([C@@H](C(=O)N[C@H]([C@@H]2C[C@H](C([C@H](O2)C[C@@H](COC)OC)(C)C)O)OC)O)OC)C",
        1.0,
        ["2-hydroxyglycine", "A5", "A8", "B1", "B2", "B3", "C1", "acetic acid", "methanol", "methylation", "methylation", "methylation", "methylation"],
    ),
    (
        "peluriside A",
        r"CC[C@@H](CO)/C=C(/C)\[C@@H]1C[C@H](C[C@@H](C([C@@]2([C@@H]([C@@H](C[C@@H](O2)C[C@H]([C@@H](C(=O)O1)O)OC)OC)O)O)(C)C)O)OC",
        1.0,
        ["A5", "B1", "B1", "B1", "B1", "B3", "B5", "C2", "D7", "acetic acid", "methylation", "methylation", "methylation"],
    ),
    (
        "penicillin G",
        r"CC1([C@@H](N2[C@H](S1)[C@@H](C2=O)NC(=O)CC3=CC=CC=C3)C(=O)O)C",
        1.0,
        ["2-phenylacetic acid", "cysteine", "valine"],
    ),
    (
        "periconiasin A",
        r"C[C@H]1[C@H]2[C@@H](NC(=O)[C@@]23[C@@H](CC(=CC[C@@H](CC3=O)O)C)C=C1C)CC(C)C",
        1.0,
        ["A1", "B1", "C1", "C1", "C2", "D2", "acetic acid", "leucine"],
    ),
    (
        "periconiasin I",
        r"C/C/1=C/C[C@@H](CC(=O)[C@]23[C@@H](C1)[C@H](C(=C([C@H]2[C@@H](NC3=O)CC(C)C)C)C)O)O",
        1.0,
        ["A1", "B1", "C1", "C1", "C2", "D2", "acetic acid", "leucine", "oxidation"],
    ),
    (
        "ratjadon",
        r"C/C=C/[C@H]1[C@H]([C@@H](C[C@H](O1)[C@@H](/C=C/C=C(\C)/C[C@@H](C)/C=C(/C)\C=C\[C@H]2CC=CC(=O)O2)O)O)C",
        1.0,
        ["B1", "B1", "B2", "C1", "C1", "C1", "C1", "C2", "C2", "D11", "D2", "acetic acid"],
    ),
    (
        "soraphen A",
        r"C[C@H]1/C=C/[C@H]([C@H](CCCC[C@H](OC(=O)[C@H]([C@@]2([C@@H]([C@H]([C@@H]([C@H]1O2)C)O)OC)O)C)C3=CC=CC=C3)OC)OC",
        1.0,
        ["A2", "B1", "B1", "B2", "B5", "C2", "D1", "D5", "benzoic acid", "methylation", "methylation", "methylation"],
    ),
    (
        "spiculoic acid A",
        r"CC[C@@H]1[C@@H]2[C@@H]([C@H](C1=O)C)C(=C[C@@]([C@@]2(CC)C(=O)O)(CC)/C=C/C3=CC=CC=C3)CC",
        1.0,
        ["2-phenylacetic acid", "A4", "C2", "C4", "C4", "C4"],
    ),
    (
        "spongidepsin",
        r"C[C@H]1CCC(CC(OC(=O)[C@@H](N(C(=O)[C@H](C1)C)C)CC2=CC=CC=C2)CCCC#C)C",
        1.0,
        ["3-butynoic acid", "B1", "D1", "D2", "D2", "D8", "methylation", "phenylalanine"],
    ),
    (
        "thailanstatin A",
        r"C[C@H]1C[C@H]([C@H](O[C@H]1C/C=C(\C)/C=C/[C@@H]2[C@H]([C@@]3(C[C@H](O2)CC(=O)O)CO3)O)C)NC(=O)/C=C\[C@H](C)OC(=O)C",
        1.0,
        ["A8", "B5", "C1", "C1", "C1", "C1", "C2", "D2", "acetic acid", "lactic acid", "oxidation", "threonine"],
    ),
    (
        "theopederin A",
        r"C[C@H]1[C@H](O[C@](CC1=C)([C@H](C(=O)N[C@@H]2[C@@H]3[C@@H]([C@H](C([C@H](O3)C[C@H]4CCCC(O4)O)(C)C)OC)OCO2)O)OC)C",
        1.0,
        ["2-hydroxyglycine", "A5", "A8", "B1", "B2", "B3", "B5", "C1", "D1", "acetic acid", "formaldehyde", "methylation", "methylation"]
    ),
    (
        "theopederin B",
        r"C[C@H]1[C@H](O[C@](CC1=C)([C@H](C(=O)N[C@@H]2[C@@H]3[C@@H]([C@H](C([C@H](O3)C[C@@H](CCCC(=O)OC)O)(C)C)OC)OCO2)O)OC)C",
        1.0,
        ["2-hydroxyglycine", "A5", "A8", "B1", "B2", "B3", "B5", "C1", "D1", "acetic acid", "formaldehyde", "methylation", "methylation", "methylation"]
    ),
    (
        "thermolide A",
        r"C[C@@H]1C[C@H]([C@@H](OC(=O)[C@H](NC(=O)C[C@H](C[C@@H]1O)O)C)[C@@H](C)C[C@H](C)[C@@H]([C@H](C)[C@@H](C[C@H](C)O)OC(=O)C)O)C",
        1.0,
        ["B1", "B1", "B1", "B2", "B2", "B2", "D2", "D2", "acetic acid", "acetic acid", "alanine"],
    ),
    (
        "thiocoraline",
        r"CN1C2CSSCC(C(=O)N(C(C(=O)SCC(C(=O)NCC1=O)NC(=O)C3=NC4=CC=CC=C4C=C3O)CSC)C)N(C(=O)CNC(=O)C(CSC(=O)C(N(C2=O)C)CSC)NC(=O)C5=NC6=CC=CC=C6C=C5O)C",
        1.0,
        ["3-hydroxyquinaldic acid", "3-hydroxyquinaldic acid", "cysteine", "cysteine", "cysteine", "cysteine", "cysteine", "cysteine", "glycine", "glycine", "methylation", "methylation", "methylation", "methylation", "methylation", "methylation"],
    ),
    (
        "zincophorin",
        r"CCC[C@@H](C)/C=C(\C)/[C@@H]([C@H](C)/C=C/CC[C@H]([C@H](C)[C@@H]([C@H](C)[C@@H]([C@H](C)[C@@H]1[C@H](CC[C@H](O1)[C@H](C)C(=O)O)C)O)O)O)O",
        1.0,
        ["B2", "B2", "B2", "B2", "B2", "C1", "C2", "C2", "D1", "D1", "D2", "propanoic acid"],
    ),
    (
        "zwittermicin A",
        r"C([C@H]([C@H](CO)N)O)[C@H]([C@H]([C@H]([C@@H](C(=O)N[C@@H](CNC(=O)N)C(=O)N)O)O)N)O",
        1.0,
        ["B1", "B12", "B5", "carbamic acid", "serine", "2,3-diaminopropionate", "amination"],
    ),
    (
        "enterobactin",
        r"C1C(C(=O)OCC(C(=O)OCC(C(=O)O1)NC(=O)C2=C(C(=CC=C2)O)O)NC(=O)C3=C(C(=CC=C3)O)O)NC(=O)C4=C(C(=CC=C4)O)O",
        1.0,
        ["2,3-dihydroxybenzoic acid", "2,3-dihydroxybenzoic acid", "2,3-dihydroxybenzoic acid", "serine", "serine", "serine"],
    ),
    (
        "curvularide C",
        r"CC[C@H](C)[C@@H](CO)NC(=O)/C=C/[C@](C)([C@H]([C@@H](C)C[C@@H](CC)O)O)OC",
        1.0,
        ["B6", "C1", "D2", "D5", "acetic acid", "isoleucinol", "methylation"],
    ),
    (
        "neopeltolide",
        r"O=C(/C=C\CCC1=COC(/C=C\CNC(OC)=O)=N1)O[C@@H]2C[C@@H](C[C@@H](C)C[C@H](OC)C[C@H](CCC)OC(C3)=O)O[C@@H]3C2",
        1.0,
        ["B1", "B1", "B1", "B1", "C1", "C1", "C1", "D1", "D1", "D8", "acetic acid", "carbonic acid", "glycine", "methylation", "methylation", "serine", "butanoic acid"],
    ),
    (
        "dihydroxydione",
        r"CCC(=O)[C@H](C)[C@H]([C@@H](C)C(=O)CC/C(=C\CC(/C(=C/C1=CSC(=N1)C)/C)O)/C)O",
        1.0,
        ["A2", "B1", "B2", "C2", "C2", "D1", "acetic acid", "cysteine", "propanoic acid"],
    ),
    (
        "amamistatin B",
        r"CCCCCCC[C@@H](C(C)(C)C(=O)N[C@H]1CCCCN(C1=O)O)OC(=O)[C@@H](CCCCN(C=O)O)NC(=O)C2=C(OC(=N2)C3=CC=CC=C3O)C",
        1.0,
        ["3-hydroxy-2,2-dimethyldecanoic acid", "B3", "D1", "D1", "D1", "N6-formyl-N6-hydroxylysine", "N6-hydroxylysine", "acetic acid", "butanoic acid", "hexanoic acid", "octanoic acid", "salicylic acid", "threonine"]
    ),
    (
        "nocardichelin B",
        r"CCCCCCCCCCC/C=C\C(=O)N(CCCCCNC(=O)CCC(=O)N(CCCCCNC(=O)[C@@H]1COC(=N1)C2=CC=CC=C2O)O)O",
        1.0,
        ["C1", "D1", "D1", "D1", "D1", "D1", "N-(5-aminopentyl)hydroxylamine", "N-(5-aminopentyl)hydroxylamine", "acetic acid", "butanedioic acid", "butanoic acid", "decanoic acid", "dodecanoic acid", "hexanoic acid", "octanoic acid", "salicylic acid", "serine", "tetradec-2-enoic acid"]
    ),
    (
        "borophycin",
        r"[B-]123O[C@]45O[C@H](C(C(=O)C[C@H](CC/C=C\C[C@@H](OC(=O)C(O1)[C@]6(O2)O[C@H](C(C(=O)C[C@H](CC/C=C\C[C@@H](OC(=O)C4O3)CC)O)(C)C)CC[C@H]6C)CC)O)(C)C)CC[C@H]5C",
        1.0,
        ["A3", "A3", "A5", "A5", "B1", "B1", "B1", "B1", "B1", "B1", "C1", "C1", "D1", "D1", "D2", "D2", "boronation", "propanoic acid", "propanoic acid"],
    ),
    (
        "aplasmomycin C",
        r"[B-]123O[C@]45O[C@H](C([C@H](C/C=C/[C@H]6O[C@@H]([C@H](C6)OC(=O)[C@@H](O1)[C@]7(O2)O[C@H](C([C@H](C/C=C/[C@H]8O[C@@H]([C@H](C8)OC(=O)[C@H]4O3)C)OC(=O)C)(C)C)CC[C@H]7C)C)OC(=O)C)(C)C)CC[C@H]5C",
        1.0,
        ["A5", "A5", "B1", "B1", "B1", "B1", "B3", "B3", "C1", "C1", "C1", "C1", "D2", "D2", "acetic acid", "acetic acid", "boronation", "lactic acid", "lactic acid"],
    ),
    (
        "chaetosemin G",
        r"Cc1c(O)c(C(O[C@H](C2)C)=O)c2c(Cl)c1O",
        1.0,
        ["A1", "A2", "B1", "C1", "acetic acid", "chlorination"],
    ),
    (
        "actinoquinoline B",
        r"CC(C)CC(=O)NC[C@H]1CC[C@@H]([C@@H](O1)O)NC(=O)C2=NC3=CC=CC=C3C=C2O",
        1.0,
        ["3-hydroxyquinaldic acid", "3-methylbutanoic acid", "5-hydroxylysine", "D8", "acetic acid"]
    ),
    (
        "NRP with intramolecular threonine-cysteine bridges",
        r"CC[C@H](C)[C@@H]1NC(=O)[C@H]2NC(=O)[C@H](Cc3ccc(O)cc3)NC(=O)CNC(=O)[C@H](Cc3c[nH]c4ccccc34)NC(=O)[C@H](Cc3c[nH]cn3)NC(=O)[C@H](CS[C@@H]2C)NC(=O)[C@H](C(C)C)NC(=O)[C@H]2NC(=O)[C@@H](NC(=O)[C@@H](NC(=O)CN)[C@@H](C)O)CS[C@H](C)[C@@H](C(=O)N[C@@H](CCCNC(=N)N)C(=O)N[C@@H](CO)C(=O)O)NC(=O)[C@H](CO)NC(=O)[C@H](Cc3ccccc3)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CS[C@@H]2C)NC(=O)[C@H](CC(C)C)NC1=O",
        1.0,
        ["arginine", "aspartic acid", "aspartic acid", "cysteine", "cysteine", "cysteine", "glycine", "glycine", "histidine", "isoleucine", "leucine", "phenylalanine", "serine", "serine", "threonine", "threonine", "threonine", "threonine", "tryptophan", "tyrosine", "valine"],
    ),
]
