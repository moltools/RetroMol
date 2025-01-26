
from typing import Optional
import json
import logging
import os
import typing as ty

from rdkit import Chem

from retromol.chem import Reaction, Reactant, get_default_linearization_rules, get_default_sequencing_rules
from retromol.matching import Motif, match_mol_greedily, greedy_max_set_cover
from retromol.react import preprocess_mol, sequence_mol
from retromol.forward import forward_generation


def run_retromol(
    name: str, 
    smiles: str, 
    user_supplied_linearization_rules: Optional[list] = None,
    user_supplied_sequencing_rules: Optional[list] = None,
    user_supplied_motfs: Optional[list] = None,
    use_default_linearization_rules: bool = True,
    use_default_sequencing_rules: bool = True,
    use_default_motifs: bool = True,
    out_folder: Optional[str] = None, 
    logger: Optional[logging.Logger] = None
) -> ty.Union[float, dict]:
    """Parse a SMILES string.
    """
    # parse the molecule
    reactant = Reactant(name, smiles)
    if logger: logger.debug(f"parsing molecule {reactant.name} with SMILES {Chem.MolToSmiles(reactant.mol)}")

    # parse user supplied rules
    if user_supplied_linearization_rules:
        if logger: logger.debug("loading user supplied linearization rules")
        user_supplied_linearization_rules = [Reaction(rule["name"], rule["reaction_smarts"]) for rule in user_supplied_linearization_rules]
        if logger: logger.debug(f"loaded {len(user_supplied_linearization_rules)} linearization reaction rules")

    if user_supplied_sequencing_rules:
        if logger: logger.debug("loading user supplied sequencing rules")
        user_supplied_sequencing_rules = [Reaction(rule["name"], rule["reaction_smarts"]) for rule in user_supplied_sequencing_rules]
        if logger: logger.debug(f"loaded {len(user_supplied_sequencing_rules)} sequencing reaction rules")

    # parse user supplied motifs
    if user_supplied_motfs:
        if logger: logger.debug("loading user supplied motifs")
        user_supplied_motfs = [Motif(motif["name"], motif["smiles"]) for motif in user_supplied_motfs]
        if logger: logger.debug(f"loaded {len(user_supplied_motfs)} motifs")
    
    # retrieve reaction rules
    if logger: logger.debug("loading default reaction rules")
    if use_default_linearization_rules:
        linearization_rules = get_default_linearization_rules()
        if logger: logger.debug(f"loaded {len(linearization_rules)} linearization reaction rules")
        linearization_rules += user_supplied_linearization_rules if user_supplied_linearization_rules else []
    else:
        linearization_rules = user_supplied_linearization_rules if user_supplied_linearization_rules else []
    
    if use_default_sequencing_rules:
        sequencing_rules = get_default_sequencing_rules()
        if logger: logger.debug(f"loaded {len(sequencing_rules)} sequencing reaction rules")
        sequencing_rules += user_supplied_sequencing_rules if user_supplied_sequencing_rules else []
    else:
        sequencing_rules = user_supplied_sequencing_rules if user_supplied_sequencing_rules else []

    # preprocess the molecule
    if logger: logger.debug(f"preprocessing molecule {name} with linearization rules")
    encoding_to_mol, encoding_to_rxn, reaction_graph, applied_reactions = preprocess_mol(reactant, linearization_rules, logger)
    leaf_nodes = [parent for parent, rxns in reaction_graph.items() if not rxns]
    if logger: logger.debug(f"found {len(leaf_nodes)} leaf nodes")

    # find identities for nodes
    encoding_to_identity = {}
    for encoding, mol in encoding_to_mol.items():
        mol_identity = match_mol_greedily(mol, user_supplied_motfs, use_default_motifs, logger)
        if mol_identity is not None:
            encoding_to_identity[encoding] = mol_identity
            continue

    # find best match for identified nodes
    matched_nodes =  [n for n, identity in encoding_to_identity.items() if identity is not None]
    identified_picked_nodes = greedy_max_set_cover(encoding_to_mol, matched_nodes)
    if logger: logger.debug(f"identified picked nodes: {identified_picked_nodes}")

    # gather all tags of the identified nodes
    matched_tags =  set()
    for node in identified_picked_nodes:
        node_tags = [atom.GetIsotope() for atom in encoding_to_mol[node].GetAtoms() if atom.GetIsotope() != 0]
        if logger: logger.debug(f"{node}, {Chem.MolToSmiles(encoding_to_mol[node])}, {node_tags}")
        matched_tags.update(node_tags)
        continue

    # identified_picked_nodes + unidentified_picked_nodes should be the set of all nodes from the input molecule
    all_tags = set(reactant.get_tags())

    # filter out any tags from all_tags that are not in matched_tags
    unmatched_tags = all_tags - matched_tags

    # find all nodes from the reaction tree that only contain unmatched tags, and could be sequenced
    unmatched_sequenceable_nodes = []
    encoding_to_motif_codes = {}
    for node, node_mol in encoding_to_mol.items():

        # only consider leaf nodes, as these should be the sequenceable nodes
        if node not in leaf_nodes:
            continue

        motif_codes = sequence_mol(node_mol, sequencing_rules, logger)
        encoding_to_motif_codes[node] = motif_codes
        if len(motif_codes) != 0:
            tags = [atom.GetIsotope() for atom in node_mol.GetAtoms() if atom.GetIsotope() != 0]
            if all(tag in unmatched_tags for tag in tags):
                unmatched_sequenceable_nodes.append(node)
    
    # find best match for unidentified but sequenceable nodes
    unidentified_sequenceable_picked_nodes = greedy_max_set_cover(encoding_to_mol, unmatched_sequenceable_nodes)

    # gather all tags from unidentified but sequenceable nodes
    tags_unidentified_nodes = set()
    for node in unidentified_sequenceable_picked_nodes:
        node_tags = [atom.GetIsotope() for atom in encoding_to_mol[node].GetAtoms() if atom.GetIsotope() != 0]
        tags_unidentified_nodes.update(node_tags)
    unmatched_tags = unmatched_tags - tags_unidentified_nodes

    # find all nodes from the reaction that only contain new unmatched tags, and cannot be sequenced
    unmatched_unsequenceable_nodes = []
    for node, node_mol in encoding_to_mol.items():
        tags = [atom.GetIsotope() for atom in node_mol.GetAtoms() if atom.GetIsotope() != 0]
        if all(tag in unmatched_tags for tag in tags):
            if node not in unmatched_sequenceable_nodes:
                unmatched_unsequenceable_nodes.append(node)

    # find best match for unidentified and unsequenceable nodes
    if unmatched_unsequenceable_nodes:
        unidentified_unsequenceable_picked_nodes = greedy_max_set_cover(encoding_to_mol, unmatched_unsequenceable_nodes)
    else:
        unidentified_unsequenceable_picked_nodes = []

    # collate all tags in identified picked nodes and separately in unidentified picked nodes
    tags_identified_nodes = set()
    for node in identified_picked_nodes:
        node_tags = [atom.GetIsotope() for atom in encoding_to_mol[node].GetAtoms() if atom.GetIsotope() != 0]
        tags_identified_nodes.update(node_tags)
    if logger:
        logger.debug(f"picked identified nodes: {identified_picked_nodes}")
        logger.debug(f"tags identified nodes: {tags_identified_nodes}")
        for node in identified_picked_nodes:
            logger.debug(f"{node} - {Chem.MolToSmiles(encoding_to_mol[node])} - {encoding_to_identity[node]}")
    
    tags_unidentified_sequenceable_nodes = set()
    for node in unidentified_sequenceable_picked_nodes:
        node_tags = [atom.GetIsotope() for atom in encoding_to_mol[node].GetAtoms() if atom.GetIsotope() != 0]
        tags_unidentified_sequenceable_nodes.update(node_tags)
    if logger:
        logger.debug(f"picked unidentified sequenceable nodes: {unidentified_sequenceable_picked_nodes}")
        logger.debug(f"tags unidentified sequenceable nodes: {tags_unidentified_sequenceable_nodes}")
        for node in unidentified_sequenceable_picked_nodes:
            logger.debug(f"{node} - {Chem.MolToSmiles(encoding_to_mol[node])}")

    tags_unidentified_unsequenceable_nodes = set()
    for node in unidentified_unsequenceable_picked_nodes:
        node_tags = [atom.GetIsotope() for atom in encoding_to_mol[node].GetAtoms() if atom.GetIsotope() != 0]
        tags_unidentified_unsequenceable_nodes.update(node_tags)
    if logger:
        logger.debug(f"picked unidentified unsequenceable nodes: {unidentified_unsequenceable_picked_nodes}")
        logger.debug(f"tags unidentified unsequenceable nodes: {tags_unidentified_unsequenceable_nodes}")
        for node in unidentified_unsequenceable_picked_nodes:
            logger.debug(f"{node} - {Chem.MolToSmiles(encoding_to_mol[node])}")

    # sanity checks based on tags, none of the tags sets can overlap
    assert len(tags_identified_nodes & tags_unidentified_sequenceable_nodes) == 0, "identified and unidentified sequenceable tags overlap"
    assert len(tags_identified_nodes & tags_unidentified_unsequenceable_nodes) == 0, "identified and unidentified unsequenceable tags overlap"

    # sanity checks based on tags, all sets combined should have same tags as all_tags
    assert tags_identified_nodes | tags_unidentified_sequenceable_nodes | tags_unidentified_unsequenceable_nodes == all_tags, "picked nodes have different tags than all tags"

    # find best motif coverage for the unidentified sequenceable nodes
    tags_identified_motifs = set()
    for node in unidentified_sequenceable_picked_nodes:
        if logger: logger.debug(f"calculating motif coverage for sequenceable node {node}")
        motif_codes = encoding_to_motif_codes[node]
        best_coverage_tags = set()
        for motif_code in motif_codes:
            tags_identified_motifs_motif_code = set()
            for motif in motif_code:
                motif_identity = match_mol_greedily(motif, user_supplied_motfs, use_default_motifs, logger)
                if motif_identity is not None:
                    motif_tags = [atom.GetIsotope() for atom in motif.GetAtoms() if atom.GetIsotope() != 0]
                    tags_identified_motifs_motif_code.update(motif_tags)
            if len(tags_identified_motifs_motif_code) > len(best_coverage_tags):
                best_coverage_tags = tags_identified_motifs_motif_code
        tags_identified_motifs.update(best_coverage_tags)

    # calculate coverage score, round to 3 decimal places
    coverage_score = (len(tags_identified_nodes) + len(tags_identified_motifs)) / len(all_tags) * 100
    coverage_score = round(coverage_score, 2)

    # collate all nodes that have been picked
    picked_nodes = identified_picked_nodes + unidentified_sequenceable_picked_nodes + unidentified_unsequenceable_picked_nodes

    # make sure picked sequenceable nodes are in identity mapping as "_backbone"
    for node in unidentified_sequenceable_picked_nodes:
        if node not in encoding_to_identity:
            encoding_to_identity[node] = "_backbone"

    # create output json file
    out_data = dict(
        name=name,
        smiles=Chem.MolToSmiles(reactant.mol),
        coverage_score=coverage_score,
        applied_preprocessin_rules=list(applied_reactions),
        encoding_to_smiles={encoding: Chem.MolToSmiles(encoding_to_mol[encoding]) for encoding in picked_nodes},
        encoding_to_identity={encoding: encoding_to_identity.get(encoding, None) for encoding in picked_nodes},
        encoding_to_motif_codes={
            encoding: {
                i: [
                    {
                        "identity": match_mol_greedily(
                            motif, 
                            user_supplied_motfs,
                            use_default_motifs,
                            logger
                        ),  # TODO: now identifying twice, once here and once in the loop above
                        "smiles": Chem.MolToSmiles(motif)
                    }
                    for motif in motif_code
                ]
                for i, motif_code in enumerate(encoding_to_motif_codes.get(encoding, [])) 
            }
            for encoding in unidentified_sequenceable_picked_nodes
        }   
    )

    # add edge list per identified node
    encoding_to_atom_list = {}
    for encoding, smiles in out_data["encoding_to_smiles"].items():
        mol = Chem.MolFromSmiles(smiles)
        atom_list = []
        for atom in mol.GetAtoms():
            atom_list.append(atom.GetIsotope())
        encoding_to_atom_list[encoding] = atom_list
    out_data["encoding_to_atom_list"] = encoding_to_atom_list

    color_pallette = [
        "#e69f00",
        "#56b4e9",
        "#039e73",
        "#f0e442",
        "#0072b2",
        "#d55f00",
        "#cc79a7",
        "#ceccca",
    ]

    polyketide_smarts = Chem.MolFromSmarts("[OH]SC~CC(=O)[OH]")
    cooh_smarts = Chem.MolFromSmarts("C(=O)[OH]")

    # get encoding to highlights 
    encoding_to_highlights = {}
    encoding_to_primary_sequence = {}
    for encoding, identity in out_data["encoding_to_identity"].items():
        if identity == "_backbone":
            motif_codes = out_data["encoding_to_motif_codes"][encoding].values()
            # every item in motif_code has 'identity', pick the one that has most items with identity not None
            best_motif_code = max(motif_codes, key=lambda x: sum(1 for motif in x if motif["identity"] is not None))
            highlights = []
            primary_sequence = []
            for i, motif in enumerate(best_motif_code):
                motif_highlights = []
                color = color_pallette[i % len(color_pallette)]
                primary_sequence.append({
                    "name": motif["identity"] if motif["identity"] is not None else "XX",
                    "smiles": motif["smiles"],
                    "color": color
                })
                motif_smiles = motif["smiles"]
                mol = Chem.MolFromSmiles(motif_smiles)

                # make sure carboxyl group is not included in highlight for output
                # only for last motif (so first in this sequence):
                cooh_inds = []
                # check if mol is polyketide motif:
                if mol.HasSubstructMatch(polyketide_smarts):
                    # see which atoms match the cooh_smarts:
                    matches = mol.GetSubstructMatches(cooh_smarts)
                    if matches:
                        for ind in matches[0]:
                            cooh_inds.append(ind)   

                for atom in mol.GetAtoms():
                    if atom.GetIdx() in cooh_inds:
                        continue
                    tag = atom.GetIsotope()
                    if tag > 0:
                        motif_highlights.append(tag)
                highlights.append(motif_highlights)
            encoding_to_highlights[encoding] = highlights[::-1]
            encoding_to_primary_sequence[encoding] = primary_sequence[::-1]
    out_data["encoding_to_highlights"] = encoding_to_highlights
    out_data["encoding_to_primary_sequence"] = encoding_to_primary_sequence

    for encoding, motifs in out_data["encoding_to_primary_sequence"].items():
        motif_smiles = []
        for motif in motifs:
            motif_smiles.append(motif["smiles"])
        try:
            linear_smiles = forward_generation(motif_smiles)
        except:
            continue
        out_data["encoding_to_smiles"][encoding] = linear_smiles
        

    if out_folder is None:
        return out_data

    out_file = os.path.join(out_folder, f"out.json")
    with open(out_file, 'w') as f:
        json.dump(out_data, f, indent=4)

    # TODO: generate output files
    # 1) data structure that shows the monomers that the reactant molecule is composed of
    # 2) data structure that allows reconstruction of the reactant molecule

    return coverage_score