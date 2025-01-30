from retromol.data.motifs import MOTIFS
from retromol.data.rules import SEQUENCING_RULES
from rdkit import Chem
from rdkit.Chem.rdChemReactions import ReactionFromSmarts
import typing as ty

def get_motif(name, as_mol = True):
    for motif in MOTIFS:
        if motif["name"] == name:
            if as_mol:
                return Chem.MolFromSmiles(motif["smiles"])
            return motif
        
def get_rules():
    # backward_rules = [r["reaction_smarts"] for r in SEQUENCING_RULES if r["props"]["include_in_forward"] == True]
    backward_rules = []
    for rule in SEQUENCING_RULES:
        if props := rule.get("props"):
            if props.get("include_in_forward"):
                backward_rules.append(rule["reaction_smarts"])

    forward_rules = []
    for rule in backward_rules:
        reactants, products = rule.split(">>")
        forward_rules.append(f"{products}>>{reactants}")
    forward_rules = [ReactionFromSmarts(rule) for rule in forward_rules]
    return forward_rules


def cap_forward_generated(mol):
    pattern = "[C:1][S][OH]>>[C:1][OH]"
    reaction = ReactionFromSmarts(pattern)
    result = reaction.RunReactants((mol,))
    if result:
        return result[0][0]
    return mol


def forward_generation(motifs: ty.List[str], as_name = False):
    if as_name:
        seq = []
        for motif in motifs:
            mol = get_motif(motif)
            if mol is None:
                assert False, f"Motif {motif} not found during forward generation"
            seq.append(mol)
    else:
        seq = []
        for motif in motifs:
            seq.append(Chem.MolFromSmiles(motif))
    rules = get_rules()
    # assemble seq into mol
    mol = seq[0]

    for motif in seq[1:]:
        # check for every rule if we can combine motif with mol
        results = []
        for rule in rules:
            result = rule.RunReactants((mol, motif))
            if result:
                # print(result)
                product = result[0][0]
                Chem.SanitizeMol(product)
                # print(product)
                results.append(product)
        assert len(results) == 1, f"More than one rule can be applied to {Chem.MolToSmiles(mol)} + {Chem.MolToSmiles(motif)}"
        mol = results[0]
    # mol = cap_forward_generated(mol)    
    return Chem.MolToSmiles(mol)