# -*- coding: utf-8 -*-

"""Compound inputs for assembly-graph integration tests.

Required assembly graphs, including coverage, live in
integration_retromol_assemblies.yml. Cases contain (name, SMILES).
"""

CASES: list[tuple[str, str]] = [
    (
        "10-deoxymethynolide",
        r"CC[C@@H]1[C@@H](/C=C/C(=O)[C@@H](C[C@@H]([C@@H]([C@H](C(=O)O1)C)O)C)C)C",
    ),
    (
        "13-deoxytedanolide",
        r"C/C=C\[C@H](C)[C@@H]1[C@](O1)(C)[C@H]([C@H]2COC(=O)[C@@H]([C@H]([C@@H](C(=O)[C@@H]([C@H](/C(=C/[C@@H](C(=O)CC[C@H](C2=O)C)C)/C)O)C)C)OC)O)O",
    ),
    (
        "2-deoxystreptamine",
        r"C1C(C(C(C(C1N)O)O)O)N",
    ),
    (
        "6-deoxyerythronolide B",
        r"CC[C@@H]1[C@@H]([C@@H]([C@H](C(=O)[C@@H](C[C@@H]([C@@H]([C@H]([C@@H]([C@H](C(=O)O1)C)O)C)O)C)C)C)O)C",
    ),
    (
        "AF-toxin",
        r"CCC(C)C(C(=O)OC(/C=C/C=C/C=C/C(=O)O)C1(CO1)C)OC(=O)C(C(C)(C)O)O",
    ),
    (
        "abyssomicin C",
        r"C[C@@H]1C[C@]23OC(=O)C4=C2OC1[C@H](O)C3\C=C\C(=O)[C@@H](C)C[C@@H](C)C4=O",
    ),
    (
        "atrop-abyssocymin C",
        r"CC1CC23OC(=O)C4=C2OC1C(O)C3\C=C/C(=O)C(C)CC(C)C4=O",
    ),
    (
        "aculeximycin",
        r"CCCC(O[C@H]1C[C@](C)(N)[C@H](O)[C@H](C)O1)C(C)C(O)C(CC)\C=C\C(O)C(C)C1C\C=C(C)\C(O)C(C)C(CC(O)C(C)C(O)CC2CC(O)C(O)C(O)(CC(O[C@@H]3O[C@H](C)[C@@H](O)[C@H](O[C@H]4C[C@@H](N)[C@H](O)[C@@H](C)O4)[C@H]3O[C@@H]3O[C@H](C)[C@@H](O)[C@H](O)[C@H]3O)C(C)CCC(O)CC(O)C\C=C(CC)\C(=O)O1)O2)O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]1O",
    ),
    (
        "acutiphycin",
        r"CCCCC[C@@H]1C/C=C(\[C@H](C(C(=O)[C@H](/C=C(\[C@@H]2C[C@@H](C[C@@](O2)(CC(=O)O1)O)O)/C)C)(C)C)O)/C",
    ),
    (
        "aflatoxin G1",
        r"COC1=C2C3=C(C(=O)OCC3)C(=O)OC2=C4[C@H]5C=CO[C@H]5OC4=C1",
    ),
    (
        "alternapyrone",
        r"CCC(C)/C=C(\C)/CCCC(C)/C=C(\C)/C=C(\C)/CC(C)C1=C(C(=C(C(=O)O1)C)O)C",
    ),
    (
        "amicoumacin",
        r"CC(C)C[C@@H]([C@@H]1CC2=C(C(=CC=C2)O)C(=O)O1)NC(=O)[C@H]([C@H]([C@H](CC(=O)N)N)O)O",
    ),
    (
        "amphidinolide J",
        r"CCC/C=C/[C@@H](C)[C@H]1C(/C=C\C([C@H](C=CCCC(=C)[C@H](CC(=O)O1)C)O)C)O",
    ),
    (
        "amphidinolide P",
        r"C[C@@H]1C(=C)C[C@H]2[C@H]3[C@@H](O3)CC(=C)/C=C/[C@H](OC(=O)C[C@@]1(O2)O)[C@H](C)C(=C)C",
    ),
    (
        "ansamitocin P-3 ",
        r"C[C@@H]1[C@@H]2C[C@]([C@@H](/C=C/C=C(/CC3=CC(=C(C(=C3)OC)Cl)N(C(=O)C[C@@H]([C@]4([C@H]1O4)C)OC(=O)C(C)C)C)\C)OC)(NC(=O)O2)O",
    ),
    (
        "anthracimycin",
        r"C[C@@H]1/C=C\C=C\[C@H](OC(=O)[C@@H](C(=O)/C=C(/[C@H]2[C@@H]1C=C[C@@H]3[C@@H]2CC=C(C3)C)\O)C)C",
    ),
    (
        "apoptolidin",
        r"COC[C@@H](C[C@H]1O[C@@](O)([C@H](O)[C@@H]2C[C@H](OC)[C@@H](O)CC\C=C(/C)\C=C\[C@@H](O[C@@H]3O[C@@H](C)[C@H](OC)[C@@H](O)[C@@H]3O)[C@H](C)\C=C(/C)\C=C(/C)\C=C(C)\C(=O)O2)[C@H](C)[C@@H](O)[C@H]1C)O[C@H]1C[C@](C)(O)[C@@H](O[C@H]2C[C@@H](OC)[C@H](O)[C@@H](C)O2)[C@H](C)O1",
    ),
    (
        "avilamycin A",
        r"C[C@@H]1[C@H]([C@@H](C[C@@H](O1)O[C@@H]2[C@H](OC3(C[C@@H]2O)O[C@@H]4[C@H](O[C@H](C[C@]4(O3)C)O[C@@H]5[C@H]([C@@H](O[C@@H]([C@@H]5OC)C)O[C@@H]6[C@H](O[C@H]([C@H]([C@H]6O)OC)O[C@H]7[C@@H]([C@H]8[C@H](CO7)O[C@@]9(O8)C1C([C@@]([C@H](O9)C)(C(=O)C)O)OCO1)OC(=O)C(C)C)COC)O)C)C)O)OC(=O)C1=C(C(=C(C(=C1OC)Cl)O)Cl)C",
    ),
    (
        "avilamycin C",
        r"C[C@@H]1[C@H]([C@@H](C[C@@H](O1)O[C@@H]2[C@H](OC3(C[C@H]2O)O[C@@H]4[C@H](O[C@H](C[C@]4(O3)C)O[C@@H]5[C@H]([C@@H](O[C@@H]([C@@H]5OC)C)O[C@@H]6[C@H](O[C@H]([C@H]([C@H]6O)OC)O[C@H]7[C@@H]([C@H]8[C@H](CO7)O[C@@]9(O8)[C@H]1[C@H]([C@@]([C@H](O9)C)(C(C)O)O)OCO1)OC(=O)C(C)C)COC)O)C)C)O)OC(=O)C1=C(C(=C(C(=C1OC)Cl)O)Cl)C",
    ),
    (
        "bitungolide F",
        r"CC[C@@H]1C=CC(=O)O[C@@H]1[C@H](C)CC[C@H](C[C@@H](/C=C/C=C/C2=CC=CC=C2)O)O",
    ),
    (
        "borrelidin",
        r"C[C@H]1C[C@H](C[C@@H]([C@H](/C(=C\C=C\C[C@H](OC(=O)C[C@@H]([C@H](C1)C)O)[C@@H]2CCC[C@H]2C(=O)O)/C#N)O)C)C",
    ),
    (
        "butirosin A",
        r"C1[C@@H]([C@H]([C@@H]([C@H]([C@@H]1NC(=O)[C@H](CCN)O)O)O[C@H]2[C@@H]([C@H]([C@H](O2)CO)O)O)O[C@@H]3[C@@H]([C@H]([C@@H]([C@H](O3)CN)O)O)N)N",
    ),
    (
        "butirosin B",
        r"C1[C@@H]([C@H]([C@@H]([C@H]([C@@H]1NC(=O)[C@H](CCN)O)O)O[C@H]2[C@@H]([C@@H]([C@H](O2)CO)O)O)O[C@@H]3[C@@H]([C@H]([C@@H]([C@H](O3)CN)O)O)N)N",
    ),
    (
        "calicheamicin",
        r"CCN[C@H]1CO[C@H](C[C@@H]1OC)O[C@@H]2[C@H]([C@@H]([C@H](OC2O[C@H]3C#C/C=C\C#C[C@]\4(CC(=O)C(=C3/C4=C\CSSSC)NC(=O)OC)O)C)NO[C@H]5C[C@@H]([C@@H]([C@H](O5)C)SC(=O)C6=C(C(=C(C(=C6OC)OC)O[C@H]7[C@@H]([C@@H]([C@H]([C@@H](O7)C)O)OC)O)I)C)O)O",
        # Default rule set is not able to parse the enediyne core, but should at least identify the glycosylation parts
    ),
    (
        "callystatin",
        r"CC[C@H](C)[C@H]([C@H](C)C(=O)[C@H](C)/C=C(\C)/C=C/C[C@@H](C)/C=C(/CC)\C=C\[C@H]1CC=CC(=O)O1)O",
    ),
    (
        "carolacton",
        r"C[C@@H]\1CCC[C@@H]([C@H](OC(=O)[C@@H]([C@@H](/C=C1)O)O)/C(=C/[C@@H](C)C(=O)[C@H](C)[C@@H](CC(=O)O)OC)/C)C",
        # This parsing is ambiguous... can parse from two sides and neither is correct/wrong
    ),
    (
        "chaetoglobosin A",
        r"C[C@H]\1C/C=C/[C@H]2[C@H]3[C@](O3)([C@H]([C@@H]4[C@@]2(C(=O)/C=C/C(=O)[C@@H](/C(=C1)/C)O)C(=O)N[C@H]4CC5=CNC6=CC=CC=C65)C)C",
    ),
    (
        "chichorine",
        r"CC1=C(O)C=C2C(CNC2=O)=C1OC",
    ),
    (
        "chlorotonil A",
        r"C[C@@H]1/C=C\C=C\[C@@H](OC(=O)[C@H](C(=O)C(C(=O)[C@@H]2[C@H]1C=C[C@H]3[C@H]2[C@@H](C=C(C3)C)C)(Cl)Cl)C)C",
    ),
    (
        "chlorothricin",
        r"C[C@@H](C(C(O)=O)=C1)C[C@@]2(C(O)=C(O3)C(O2)=O)[C@@H]1C=CCCCC[C@@H]4C=C[C@@]([C@@H](O[C@@H]5C[C@H](O)[C@@H](O[C@H]6O[C@@H](C)[C@H](O)[C@@H](OC(C7=C(C)C(Cl)=CC=C7OC)=O)C6)[C@@H](C)O5)CCC8)([H])[C@]8([H])[C@@H]4C3=O",
    ),
    (
        "coelimycin P1",
        r"C/C=C/C(=O)/C/1=C/C(=C\2/C=CCCN2)/SC[CH](C(=O)O1)NC(=O)C",
    ),
    (
        "compactin",
        r"CC[C@H](C)C(=O)O[C@H]1CCC=C2C=C[C@H](C)[C@H](CC[C@@H](O)C[C@@H](O)CC(O)=O)[C@@H]12",
    ),
    (
        "cremimycin",
        r"CCCCCCC1CC(CCCC(C/C(O)=C2C(/C=C(/C=C\C=C/C(N1)=O)C)CC(C\2=O)O[C@H]3C[C@@H]([C@@H]([C@H](O3)C)O)OC)O)=O",
    ),
    (
        "deschlorothricin",
        r"C[C@@H](C(C(O)=O)=C1)C[C@@]2(C(O)=C(O3)C(O2)=O)[C@@H]1C=CCCCC[C@@H]4C=C[C@@]([C@@H](O[C@@H]5C[C@H](O)[C@@H](O[C@H]6O[C@@H](C)[C@H](O)[C@@H](OC(C7=C(C)C=CC=C7OC)=O)C6)[C@@H](C)O5)CCC8)([H])[C@]8([H])[C@@H]4C3=O",
    ),
    (
        "daptomycin",
        r"CCCCCCCCCC(=O)N[C@@H](CC1=CNC2=CC=CC=C21)C(=O)N[C@H](CC(=O)N)C(=O)N[C@@H](CC(=O)O)C(=O)N[C@H]3[C@H](OC(=O)[C@@H](NC(=O)[C@@H](NC(=O)[C@H](NC(=O)CNC(=O)[C@@H](NC(=O)[C@H](NC(=O)[C@@H](NC(=O)[C@@H](NC(=O)CNC3=O)CCCN)CC(=O)O)C)CC(=O)O)CO)[C@H](C)CC(=O)O)CC(=O)C4=CC=CC=C4N)C",
    ),
    (
        "dictyostatin",
        r"C[C@H]1CC[C@H]([C@@H]([C@@H](OC(=O)/C=C\C=C\[C@H]([C@H](C[C@@H](/C=C\[C@@H]([C@@H]([C@H](C1)C)O)C)O)O)C)[C@@H](C)/C=C\C=C)C)O",
    ),
    (
        "discodermolide",
        r"C[C@H]1[C@@H](OC(=O)[C@@H]([C@H]1O)C)C[C@@H](/C=C\[C@H](C)[C@@H]([C@@H](C)/C=C(/C)\C[C@H](C)[C@H]([C@H](C)[C@H]([C@@H](C)/C=C\C=C)OC(=O)N)O)O)O",
    ),
    (
        "epothilone",
        r"C[C@H]1CCC[C@@H]2[C@@H](O2)C[C@H](OC(=O)C[C@H](C(C(=O)[C@@H]([C@H]1O)C)(C)C)O)/C(=C/C3=CSC(=N3)C)/C",
    ),
    (
        "erythromycin",
        r"CC[C@@H]1[C@@]([C@@H]([C@H](C(=O)[C@@H](C[C@@]([C@@H]([C@H]([C@@H]([C@H](C(=O)O1)C)O[C@H]2C[C@@]([C@H]([C@@H](O2)C)O)(C)OC)C)O[C@H]3[C@@H]([C@H](C[C@H](O3)C)N(C)C)O)(C)O)C)C)O)(C)O",
    ),
    (
        "georatusin",
        r"CC[C@@H](C)C[C@H](C)[C@@H]1[C@H](C[C@@H]([C@@H]2[C@H](C[C@@H]([C@@](O2)([C@H](C(=O)N[C@@H](C(=O)O1)CC3=CNC4=CC=CC=C43)C)O)C)C)C)C",
    ),
    (
        "gephyronic acid",
        r"C[C@@H]1[C@@H](O[C@@](C([C@H]1OC)(C)C)([C@@H](C)C[C@H](C)[C@@H]([C@@]2([C@H](O2)[C@@H](C)C=C(C)C)C)O)O)CC(=O)O",
    ),
    (
        "harzianic acid",
        r"CCC/C=C/C=C/C(=C\1/C(=O)C(N(C1=O)C)CC(C(C)C)(C(=O)O)O)/O",
    ),
    (
        "herboxidiene",
        r"C[C@H]1CC[C@@H](O[C@@H]1/C(=C/C=C/[C@@H](C)C[C@@]2([C@H](O2)[C@H](C)[C@H]([C@@H](C)O)OC)C)/C)CC(=O)O",
    ),
    (
        "hydroxystreptomycin",
        r"CN[C@H]1[C@H](O)[C@@H](O)[C@H](CO)O[C@H]1O[C@H]1[C@H](O[C@H]2[C@H](O)[C@@H](O)[C@H](NC(N)=N)[C@@H](O)[C@@H]2NC(N)=N)O[C@@H](CO)[C@]1(O)C=O",
    ),
    (
        "hymenosetin",
        r"C/C=C/[C@@H]1C(=C[C@@H]2C[C@@H](CC[C@H]2[C@]1(C)/C(=C\3/C(=O)[C@@H](NC3=O)[C@@H](C)O)/O)C)C",
    ),
    (
        "indanomycin",
        r"CC[C@H]1CC[C@@H]2[C@@H]1C=C[C@H]([C@H]2C(=O)C3=CC=CN3)/C=C/C=C(\CC)/[C@H]4[C@H](CC[C@@H](O4)[C@@H](C)C(=O)O)C",
    ),
    (
        "ircinianin",
        r"C[C@H]1CCC2[C@@H]1C3(C(C=C2C)/C=C(\C)/CCCC4=COC=C4)C(=C(C(=O)O3)C)O",
    ),
    (
        "iriomoteolide 1a",
        r"C[C@H]1C/C=C/[C@@]([C@@]2(CC(=C)C[C@@H](O2)C/C=C/[C@@H]([C@@H](/C(=C\C(=O)O[C@@H]1C[C@H](C)[C@H](C)O)/C)C)O)O)(C)O",
    ),
    (
        "iriomoteolide 3a",
        r"C[C@H]1C/C=C/[C@@]([C@@]2(CC(=C)C[C@@H](O2)C/C=C/[C@@H]([C@@H](/C(=C\C(=O)O[C@@H]1C[C@H](C)[C@H](C)O)/C)C)O)O)(C)O",
    ),
    (
        "jerangolid A",
        r"CC[C@@H]1C(=CC[C@@H](O1)/C(=C/[C@H](C)/C=C/[C@H]2CC(=C(C(=O)O2)CO)OC)/C)C",
    ),
    (
        "kirromycin",
        r"CC[C@H](C(=O)NC/C=C/C=C(\C)/[C@H]([C@@H](C)[C@H]1[C@H]([C@H]([C@H](O1)/C=C/C=C/C=C(\C)/C(=O)C2=C(C=CNC2=O)O)O)O)OC)[C@@]3([C@@H]([C@@H](C([C@@H](O3)/C=C/C=C/C)(C)C)O)O)O",
    ),
    (
        "lactimidomycin",
        r"C[C@H]1/C=C\C=C\CC/C=C/C(=O)O[C@H]1/C(=C/[C@H](C)C(=O)C[C@@H](CC2CC(=O)NC(=O)C2)O)/C",
    ),
    (
        "lankamycin",
        r"C[C@@H]1C[C@@H]([C@H]([C@@H](O1)O[C@H]2[C@H](C[C@](C(=O)[C@@H]([C@H]([C@H]([C@H](OC(=O)[C@@H]([C@H]([C@@H]2C)O[C@H]3C[C@@]([C@@H]([C@@H](O3)C)OC(=O)C)(C)OC)C)[C@@H](C)[C@H](C)O)C)OC(=O)C)C)(C)O)C)O)OC",
    ),
    (
        "latrunculin",
        r"C[C@H]1CC[C@@H]2C[C@H](C[C@@](O2)([C@@H]3CSC(=O)N3)O)OC(=O)C=C(CCC=CC=C1)C",
    ),
    (
        "leiodermatolide",
        r"CC[C@H]1[C@@H]([C@@](CC(=O)O1)(C/C=C/C=C(\C)/[C@@H]2[C@H](/C=C\C=C/[C@@H]([C@H]([C@H]([C@@H](/C(=C/CCC(=O)O2)/C)C)O)C)OC(=O)N)C)O)C",
    ),
    (
        "lovastatin",
        r"CC[C@H](C)C(=O)O[C@H]1C[C@H](C=C2[C@H]1[C@H]([C@H](C=C2)C)CC[C@@H]3C[C@H](CC(=O)O3)O)C",
    ),
    (
        "macrolactin A",
        r"CC1CCCC=CC=CC(CC(CC=CC=CC(CC=CC=CC(=O)O1)O)O)O",
    ),
    (
        "maytansine",
        r"C[C@@H]1[C@@H]2C[C@]([C@@H](/C=C/C=C(/CC3=CC(=C(C(=C3)OC)Cl)N(C(=O)C[C@@H]([C@]4([C@H]1O4)C)OC(=O)[C@H](C)N(C)C(=O)C)C)\C)OC)(NC(=O)O2)O",
    ),
    (
        "megalomycin A",
        r"CC[C@@H]1[C@@]([C@@H]([C@H](C(=O)[C@@H](C[C@@]([C@@H]([C@H]([C@@H]([C@H](C(=O)O1)C)O[C@H]2C[C@@]([C@H]([C@@H](O2)C)O)(C)O)C)O[C@H]3[C@@H]([C@H](C[C@H](O3)C)N(C)C)O)(C)O[C@H]4C[C@H]([C@H]([C@@H](O4)C)O)N(C)C)C)C)O)(C)O",
    ),
    (
        "micacocidin A",
        r"CCCCCC1=C(C(=CC=C1)[O-])C2=N[C@H](CS2)[C@@H]3N([C@@H](CS3)[C@@H](C(C)(C)C4=N[C@@](CS4)(C)C(=O)[O-])O)C",
    ),
    (
        "migrastatin",
        r"C[C@@H]1/C=C(\[C@H](OC(=O)/C=C/CC/C=C/[C@@H]([C@H]1O)OC)[C@H](C)C(=O)CCCC2CC(=O)NC(=O)C2)/C",
    ),
    (
        "narbonolide",
        r"CC[C@@H]1[C@@H](/C=C/C(=O)[C@@H](C[C@@H]([C@@H]([C@H](C(=O)[C@H](C(=O)O1)C)C)O)C)C)C",
    ),
    (
        "pederin",
        r"C[C@H]1[C@H](O[C@](CC1=C)([C@@H](C(=O)N[C@H]([C@@H]2C[C@H](C([C@H](O2)C[C@@H](COC)OC)(C)C)O)OC)O)OC)C",
    ),
    (
        "peluriside A",
        r"CC[C@@H](CO)/C=C(/C)\[C@@H]1C[C@H](C[C@@H](C([C@@]2([C@@H]([C@@H](C[C@@H](O2)C[C@H]([C@@H](C(=O)O1)O)OC)OC)O)O)(C)C)O)OC",
    ),
    (
        "penicillin G",
        r"CC1([C@@H](N2[C@H](S1)[C@@H](C2=O)NC(=O)CC3=CC=CC=C3)C(=O)O)C",
    ),
    (
        "periconiasin A",
        r"C[C@H]1[C@H]2[C@@H](NC(=O)[C@@]23[C@@H](CC(=CC[C@@H](CC3=O)O)C)C=C1C)CC(C)C",
    ),
    (
        "periconiasin I",
        r"C/C/1=C/C[C@@H](CC(=O)[C@]23[C@@H](C1)[C@H](C(=C([C@H]2[C@@H](NC3=O)CC(C)C)C)C)O)O",
    ),
    (
        "ratjadon",
        r"C/C=C/[C@H]1[C@H]([C@@H](C[C@H](O1)[C@@H](/C=C/C=C(\C)/C[C@@H](C)/C=C(/C)\C=C\[C@H]2CC=CC(=O)O2)O)O)C",
    ),
    (
        "soraphen A",
        r"C[C@H]1/C=C/[C@H]([C@H](CCCC[C@H](OC(=O)[C@H]([C@@]2([C@@H]([C@H]([C@@H]([C@H]1O2)C)O)OC)O)C)C3=CC=CC=C3)OC)OC",
    ),
    (
        "spiculoic acid A",
        r"CC[C@@H]1[C@@H]2[C@@H]([C@H](C1=O)C)C(=C[C@@]([C@@]2(CC)C(=O)O)(CC)/C=C/C3=CC=CC=C3)CC",
    ),
    (
        "spongidepsin",
        r"C[C@H]1CCC(CC(OC(=O)[C@@H](N(C(=O)[C@H](C1)C)C)CC2=CC=CC=C2)CCCC#C)C",
    ),
    (
        "thailanstatin A",
        r"C[C@H]1C[C@H]([C@H](O[C@H]1C/C=C(\C)/C=C/[C@@H]2[C@H]([C@@]3(C[C@H](O2)CC(=O)O)CO3)O)C)NC(=O)/C=C\[C@H](C)OC(=O)C",
    ),
    (
        "theopederin A",
        r"C[C@H]1[C@H](O[C@](CC1=C)([C@H](C(=O)N[C@@H]2[C@@H]3[C@@H]([C@H](C([C@H](O3)C[C@H]4CCCC(O4)O)(C)C)OC)OCO2)O)OC)C",
    ),
    (
        "theopederin B",
        r"C[C@H]1[C@H](O[C@](CC1=C)([C@H](C(=O)N[C@@H]2[C@@H]3[C@@H]([C@H](C([C@H](O3)C[C@@H](CCCC(=O)OC)O)(C)C)OC)OCO2)O)OC)C",
    ),
    (
        "thermolide A",
        r"C[C@@H]1C[C@H]([C@@H](OC(=O)[C@H](NC(=O)C[C@H](C[C@@H]1O)O)C)[C@@H](C)C[C@H](C)[C@@H]([C@H](C)[C@@H](C[C@H](C)O)OC(=O)C)O)C",
    ),
    (
        "thiocoraline",
        r"CN1C2CSSCC(C(=O)N(C(C(=O)SCC(C(=O)NCC1=O)NC(=O)C3=NC4=CC=CC=C4C=C3O)CSC)C)N(C(=O)CNC(=O)C(CSC(=O)C(N(C2=O)C)CSC)NC(=O)C5=NC6=CC=CC=C6C=C5O)C",
    ),
    (
        "zincophorin",
        r"CCC[C@@H](C)/C=C(\C)/[C@@H]([C@H](C)/C=C/CC[C@H]([C@H](C)[C@@H]([C@H](C)[C@@H]([C@H](C)[C@@H]1[C@H](CC[C@H](O1)[C@H](C)C(=O)O)C)O)O)O)O",
    ),
    (
        "zwittermicin A",
        r"C([C@H]([C@H](CO)N)O)[C@H]([C@H]([C@H]([C@@H](C(=O)N[C@@H](CNC(=O)N)C(=O)N)O)O)N)O",
    ),
    (
        "enterobactin",
        r"C1C(C(=O)OCC(C(=O)OCC(C(=O)O1)NC(=O)C2=C(C(=CC=C2)O)O)NC(=O)C3=C(C(=CC=C3)O)O)NC(=O)C4=C(C(=CC=C4)O)O",
    ),
    (
        "curvularide C",
        r"CC[C@H](C)[C@@H](CO)NC(=O)/C=C/[C@](C)([C@H]([C@@H](C)C[C@@H](CC)O)O)OC",
    ),
    (
        "neopeltolide",
        r"O=C(/C=C\CCC1=COC(/C=C\CNC(OC)=O)=N1)O[C@@H]2C[C@@H](C[C@@H](C)C[C@H](OC)C[C@H](CCC)OC(C3)=O)O[C@@H]3C2",
    ),
    (
        "dihydroxydione",
        r"CCC(=O)[C@H](C)[C@H]([C@@H](C)C(=O)CC/C(=C\CC(/C(=C/C1=CSC(=N1)C)/C)O)/C)O",
    ),
    (
        "amamistatin B",
        r"CCCCCCC[C@@H](C(C)(C)C(=O)N[C@H]1CCCCN(C1=O)O)OC(=O)[C@@H](CCCCN(C=O)O)NC(=O)C2=C(OC(=N2)C3=CC=CC=C3O)C",
    ),
    (
        "nocardichelin B",
        r"CCCCCCCCCCC/C=C\C(=O)N(CCCCCNC(=O)CCC(=O)N(CCCCCNC(=O)[C@@H]1COC(=N1)C2=CC=CC=C2O)O)O",
    ),
    (
        "borophycin",
        r"[B-]123O[C@]45O[C@H](C(C(=O)C[C@H](CC/C=C\C[C@@H](OC(=O)C(O1)[C@]6(O2)O[C@H](C(C(=O)C[C@H](CC/C=C\C[C@@H](OC(=O)C4O3)CC)O)(C)C)CC[C@H]6C)CC)O)(C)C)CC[C@H]5C",
    ),
    (
        "aplasmomycin C",
        r"[B-]123O[C@]45O[C@H](C([C@H](C/C=C/[C@H]6O[C@@H]([C@H](C6)OC(=O)[C@@H](O1)[C@]7(O2)O[C@H](C([C@H](C/C=C/[C@H]8O[C@@H]([C@H](C8)OC(=O)[C@H]4O3)C)OC(=O)C)(C)C)CC[C@H]7C)C)OC(=O)C)(C)C)CC[C@H]5C",
    ),
    (
        "chaetosemin G",
        r"Cc1c(O)c(C(O[C@H](C2)C)=O)c2c(Cl)c1O",
    ),
    (
        "actinoquinoline B",
        r"CC(C)CC(=O)NC[C@H]1CC[C@@H]([C@@H](O1)O)NC(=O)C2=NC3=CC=CC=C3C=C2O",
    ),
    (
        "NRP with intramolecular threonine-cysteine bridges",
        r"CC[C@H](C)[C@@H]1NC(=O)[C@H]2NC(=O)[C@H](Cc3ccc(O)cc3)NC(=O)CNC(=O)[C@H](Cc3c[nH]c4ccccc34)NC(=O)[C@H](Cc3c[nH]cn3)NC(=O)[C@H](CS[C@@H]2C)NC(=O)[C@H](C(C)C)NC(=O)[C@H]2NC(=O)[C@@H](NC(=O)[C@@H](NC(=O)CN)[C@@H](C)O)CS[C@H](C)[C@@H](C(=O)N[C@@H](CCCNC(=N)N)C(=O)N[C@@H](CO)C(=O)O)NC(=O)[C@H](CO)NC(=O)[C@H](Cc3ccccc3)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CS[C@@H]2C)NC(=O)[C@H](CC(C)C)NC1=O",
    ),
    (
        "vancomycin",
        r"CN[C@H](CC(C)C)C(=O)N[C@H]1C(=O)N[C@@H](CC(N)=O)C(=O)N[C@H]2C(=O)N[C@H]3C(=O)N[C@H](C(=O)N[C@H](C(=O)O)c4cc(O)cc(O)c4-c4cc3ccc4O)[C@H](O)c3ccc(c(Cl)c3)Oc3cc2cc(c3O[C@@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O[C@H]2C[C@](C)(N)[C@H](O)[C@H](C)O2)Oc2ccc(cc2Cl)[C@H]1O",
    ),
    (
        "capreomycin IA",
        r"NCCC[C@@H](N)CC(=O)NC[C@@H]1NC(=O)[C@H](CO)NC(=O)[C@@H](N)CNC(=O)[C@H]([C@H]2CCN=C(N)N2)NC(=O)/C(=C/NC(N)=O)NC1=O",
    ),
    (
        "aeruginosin 98-A",
        r"CC[C@@H](C)[C@@H](NC(=O)C(O)Cc1ccc(O)c(Cl)c1)C(=O)N1C(C(=O)NCCCCNC(=N)N)CC2CCC(OS(=O)(=O)O)CC21",
    ),
    (
        "corpeptin A",
        r"C=C(N=C(O)[C@@H](N=C(O)[C@@H](N=C(O)[C@H](CCO)N=C(O)/C(=C/C)N=C(O)[C@@H](N=C(O)[C@@H](N=C(O)[C@H](C)N=C(O)[C@H](C)N=C(O)[C@H](C)N=C(O)[C@@H]1CCCN1C(=O)/C(=C/C)N=C(O)CC(O)CCCCCCC)C(C)C)C(C)C)C(C)C)[C@H](C)CC)C(O)=N[C@@H](C)C(O)=N[C@@H](C)C(O)=N[C@@H](C)C(O)=N[C@H](C(O)=N/C(=C\C)C(O)=N[C@@H]1C(O)=N[C@@H](C)C(O)=N[C@@H](CCN)C(O)=N[C@@H](CO)C(O)=N[C@@H]([C@@H](C)CC)C(=O)O[C@H]1C)C(C)C",
    ),
    (
        "cyanopeptolin 1138",
        r"CCC(C)C1NC(=O)C(Cc2ccccc2)N(C)C(=O)C(C(C)CC)N2C(=O)C(CCC2O)NC(=O)C(CCc2ccc(O)cc2)NC(=O)C(NC(=O)C(CCc2ccc(O)cc2)NC(=O)C(COS(=O)(=O)O)OC)C(C)OC1=O",
    ),
    (
        "gramicidin",
        r"CC(C)C[C@@H](NC(=O)[C@H](C)NC(=O)CNC(=O)[C@@H](NC=O)C(C)C)C(=O)N[C@@H](C)C(=O)N[C@@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@@H](Cc1c[nH]c2ccccc12)C(=O)N[C@H](CC(C)C)C(=O)N[C@@H](Cc1c[nH]c2ccccc12)C(=O)N[C@H](CC(C)C)C(=O)N[C@@H](Cc1c[nH]c2ccccc12)C(=O)N[C@H](CC(C)C)C(=O)N[C@@H](Cc1c[nH]c2ccccc12)C(=O)NCCO)C(C)C)C(C)C)C(C)C",
    ),
    (
        "hormaomycin",
        r"C/C=C\[C@H]1C[C@H]2C(=O)O[C@H](C)[C@@H](NC(=O)[C@H](C[C@@H]3C[C@H]3[N+](=O)[O-])NC(=O)c3ccc(Cl)n3O)C(=O)N[C@@H]([C@H](C)c3ccccc3)C(=O)N[C@H](C[C@@H]3C[C@H]3[N+](=O)[O-])C(=O)N[C@@H]([C@H](C)c3ccccc3)C(=O)N[C@@H]([C@@H](C)CC)C(=O)N2C1",
    ),
    (
        "myxochelin A",
        r"O=C(NCCCC[C@@H](CO)NC(=O)c1cccc(O)c1O)c1cccc(O)c1O",
    ),
    (
        "marformycin A",
        r"CC[C@H](C)[C@H]1NC(=O)[C@@H](Cc2ccc(OC)cc2)NC(=O)[C@@H](NC(=O)[C@@H](NC=O)C(C)C)[C@@H](C)OC(=O)[C@H](C(C)C)N(C)C(=O)[C@H](CC(C)C)NC(=O)[C@@H]2CCCNN2C1=O",
    ),
    (
        "antimycin",
        r"CCCCCC[C@H]1C(=O)O[C@H](C)[C@H](NC(=O)c2cccc(NC=O)c2O)C(=O)O[C@@H](C)[C@@H]1OC(=O)CC(C)C",
    ),
    (
        "phorbactazole A",
        r"CO[C@@H]1C[C@H](C)OC(O)(Cc2nc(/C=C(\C)[C@@H]3O[C@@H]4C/C=C(\C)c5nc(co5)[C@@H]5C[C@@H](OC(N)=O)[C@H](O)[C@H](C[C@H](O)C[C@@H](O)[C@@H](O)/C=C/[C@@H]6C[C@@H](O)[C@H](C)[C@@H](C/C=C\C(=O)O[C@@H](C4)[C@H]3C)O6)O5)co2)[C@H]1C",
    ),
]
