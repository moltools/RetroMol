from copy import deepcopy
from typing import Dict, List, Optional, Tuple
import itertools
import logging

from rdkit import Chem
import networkx as nx

from retromol.chem import encode_mol, neutralize_mol
from retromol.react import Reaction, get_default_linearization_rules, get_default_sequencing_rules, get_default_motifs


def preprocess_mol(
    smiles: str, 
    reaction_rules: List[Reaction],
    logger: logging.Logger
) -> Tuple[Dict[int, Chem.Mol], Dict[int, Reaction], Dict[int, Dict[int, List[int]]]]:
    """Apply custom rules to linearize a SMILES string.
    
    :param smiles: SMILES string
    :type smiles: str
    :param reaction_rules: List of reaction rules
    :type reaction_rules: List[Reaction]
    :param logger: Logger object
    :type logger: logging.Logger
    :return: A tuple containing the following:
        - A dictionary mapping molecule encodings to RDKit molecules.
        - A dictionary mapping reaction encodings to Reaction objects.
        - A dictionary representing the reaction graph.
    """
    # parse the input SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"invalid input SMILES: {smiles}")
    
    # neutralize the molecule
    neutralize_mol(mol)

    # sort reaction rules based on number of heavy atoms in reactants
    # reactions with less heavy atoms should be applied first
    reaction_rule_sizes = []
    for reaction_rule_idx, reaction_rule in enumerate(reaction_rules):
        reactants = [reaction_rule.rxn.GetReactantTemplate(i) for i in range(reaction_rule.rxn.GetNumReactantTemplates())]
        num_heavy_atoms = sum([reactant.GetNumHeavyAtoms() for reactant in reactants])
        reaction_rule_sizes.append((reaction_rule_idx, num_heavy_atoms))
    reaction_rule_sizes = sorted(reaction_rule_sizes, key=lambda x: x[0])
    sorted_reaction_rules = [reaction_rules[idx] for idx, _ in reaction_rule_sizes]

    # prioritize 1-to-1 reaction rules over 1-to-many rules
    sorted_reaction_rules = sorted(sorted_reaction_rules, key=lambda x: x.rxn.GetNumProductTemplates())

    # store atom indices as isotope numbers, since those persist through reactions
    for atom in mol.GetAtoms():
        tag = atom.GetIdx() + 1
        atom.SetIsotope(tag)
    original_taken_tags = [atom.GetIsotope() for atom in mol.GetAtoms() if atom.GetIsotope() != 0]

    # create data structures
    encoding_to_mol = {}  # maps mol hash to Chem.Mol
    num_rxn_nodes = 0  # keeps track of number of reaction nodes in reaction graph
    encoding_to_rxn = {}  # maps rxn index to Reaction
    reaction_graph = {}  # as {reactant_mol_encoding: {reaction_encoding: [child_mol_encodings], ...}, ...}

    # create queue
    mols = [mol]

    # process queue
    while mols:
        parent = mols.pop(0)

        # encode parent molecule
        parent_encoding = encode_mol(parent)
        if parent_encoding not in encoding_to_mol:
            encoding_to_mol[parent_encoding] = deepcopy(parent)
            reaction_graph[parent_encoding] = {}

        # first loop over reactions and check which reaction templates are uncontested
        # apply these uncontested reactions in bulk and keep looping until no more uncontested reactions can be applied
        up_for_election = []
        for reaction_rule in sorted_reaction_rules:
            all_reactant_matches = []
            for template_idx in range(reaction_rule.rxn.GetNumReactantTemplates()):
                reactant_matches = parent.GetSubstructMatches(reaction_rule.rxn.GetReactantTemplate(template_idx))
                all_reactant_matches.append(reactant_matches)

            # generate all possible match sets, for when reaction template has multiple reactants
            match_sets = list(itertools.product(*all_reactant_matches))

            # collate all match sets into single sets of atoms
            match_sets = [set(itertools.chain(*match_set)) for match_set in match_sets]

            # add match sets to up_for_election
            for match_set in match_sets:
                up_for_election.append((reaction_rule, match_set))

        # check which reactions with matched templates are uncontested and which are contested
        uncontested = []
        contested = []
        for i, (reaction_rule, match_set) in enumerate(up_for_election):
            if not any(match_set.intersection(match_set_other) for j, (_, match_set_other) in enumerate(up_for_election) if i != j):
                uncontested.append((reaction_rule, match_set))
            else:
                contested.append((reaction_rule, match_set))

        # if neither uncontested or contested, start over with new parent
        if not uncontested and not contested:
            continue

        # apply uncontested reactions
        if uncontested:
            # we make sure all atoms, even the ones not from original reactant, have a
            # unique isotope number, so we can track them through reactions
            temp_taken_tags = [atom.GetIsotope() for atom in parent.GetAtoms() if atom.GetIsotope() != 0]
            for atom in parent.GetAtoms():
                if atom.GetIsotope() == 0:
                    tag = 1
                    while tag in temp_taken_tags:
                        tag += 1
                    atom.SetIsotope(tag)
                    temp_taken_tags.append(tag)
            
            # we don't have to worry about tagging new atoms, as all of these reactions
            # are now uncontested, meaning they don't share any atoms with other reactions
            # do make sure every atom in input is tagged
            tagged_atoms = set([atom.GetIsotope() for atom in parent.GetAtoms() if atom.GetIsotope() != 0])
            if len(tagged_atoms) != len(parent.GetAtoms()):
                raise ValueError("not all atoms have a unique tag before applying uncontested reactions")

            # map tags to atomic nums so we can reset later after using mask
            idx_to_tag = {atom.GetIdx(): atom.GetIsotope() for atom in parent.GetAtoms()}

            # all uncontested reactions become one noe in the reaction_graph teogether
            products = []

            # apply all uncontested reactions 
            for reaction_rule, match_set in uncontested:
                
                # we use the isotope numbers to mask all atoms except the ones in the match set
                mask = set([idx_to_tag[idx] for idx in match_set])
                
                # we use the original parent if there are no products, otherwise we have to use the correct
                # product from all products -- the one that contains our mask
                if len(products) != 0:
                    for product in products:
                        product_tags = set([atom.GetIsotope() for atom in product.GetAtoms() if atom.GetIsotope() != 0])
                        if mask.issubset(product_tags):
                            parent = product
                            products = [p for p in products if p != product]
                            break
                
                # make sure all atoms in parent have a unique tag
                temp_taken_tags_uncontested = [atom.GetIsotope() for atom in parent.GetAtoms() if atom.GetIsotope() != 0]
                for atom in parent.GetAtoms():
                    if atom.GetIsotope() == 0:
                        tag = 1
                        while tag in temp_taken_tags_uncontested:
                            tag += 1
                        atom.SetIsotope(tag)
                        temp_taken_tags_uncontested.append(tag)

                # apply reaction rule
                results = reaction_rule(parent, mask=mask, logger=logger)

                # we except a single result
                if len(results) == 0:
                    raise ValueError(f"no products from uncontested reaction {reaction_rule.name}")
                if len(results) > 1:
                    raise ValueError(f"more than one product from uncontested reaction {reaction_rule.name}")
                result = results[0]

                # reset atom tags in products for atoms not in original reactant
                for product in result:
                    for atom in product.GetAtoms():
                        if atom.GetIsotope() not in original_taken_tags and atom.GetIsotope() != 0:
                            atom.SetIsotope(0)
                    products.append(product)
            
            # all products are now products of our combined reaction node that
            # contains all uncontested reactions
            num_rxn_nodes += 1
            reaction_graph[parent_encoding][num_rxn_nodes] = set()
            encoding_to_rxn[num_rxn_nodes] = "uncontested"
            for product in products:
                product_encoding = encode_mol(product)
                mols.append(product)
                reaction_graph[parent_encoding][num_rxn_nodes].add(product_encoding)

            # restart loop with new parent
            continue

        # check if a 1-to-1 reaction is contested, if so, apply these first before doing anything else
        temp_reaction_rules = []
        if contested:
            for reaction_rule, match_set in contested:
                if reaction_rule.rxn.GetNumReactantTemplates() == 1 and reaction_rule.rxn.GetNumProductTemplates() == 1:
                    temp_reaction_rules.append(reaction_rule)
        if not temp_reaction_rules:
            temp_reaction_rules= sorted_reaction_rules

        # exhaustive reaction_rule application for all contested reactions
        for reaction_rule in temp_reaction_rules:
            results = reaction_rule(parent, logger=logger)
            for result in results:

                # encode reaction node
                num_rxn_nodes += 1
                if num_rxn_nodes not in encoding_to_rxn:
                    encoding_to_rxn[num_rxn_nodes] = reaction_rule.name
                else:
                    raise ValueError(f"reaction node {num_rxn_nodes} already exists for reaction {reaction_rule.name}")
                
                # encode child molecules
                reaction_graph[parent_encoding][num_rxn_nodes] = set()

                for child in result:
                    child_encoding = encode_mol(child)
                    if child_encoding not in encoding_to_mol:
                        encoding_to_mol[child_encoding] = deepcopy(child)
                        reaction_graph[child_encoding] = {}
                        mols.append(child)
                    reaction_graph[parent_encoding][num_rxn_nodes].add(child_encoding)

    # return data structures
    return encoding_to_mol, encoding_to_rxn, reaction_graph


def match_mol(mol: Chem.Mol, logger: Optional[logging.Logger] = None) -> Optional[str]:
    """Match a molecule to a motif.
    
    :param mol: RDKit molecule
    :type mol: Chem.Mol
    :param logger: Logger object
    :type logger: Optional[logging.Logger]
    :return: The name of the motif if a match is found, otherwise None
    :rtype: Optional[str]
    """
    motifs = get_default_motifs()
    for motif in motifs:
        if motif.is_match(mol):
            if logger: logger.debug(f"matched motif {motif.name} to molecule {Chem.MolToSmiles(mol)}")
            return motif.name
    
    # if no motif matched, return None
    if logger: logger.debug(f"no motif matched to molecule {Chem.MolToSmiles(mol)}")
    return None


def sequence_mol(
    mol_encoding: int,
    encoding_to_mol: Dict[int, Chem.Mol],
    encoding_to_rxn: Dict[int, Reaction],
    reaction_graph: Dict[int, Dict[int, List[int]]],
    sequencing_rules: List[Reaction],
    logger: Optional[logging.Logger] = None
) -> Tuple[Dict[int, Chem.Mol], Dict[int, Reaction], Dict[int, Dict[int, List[int]]]]:
    """Apply custom rules to sequence a molecule into motifs.
    
    :param mol_encoding: Encoding of the molecule
    :type mol_encoding: int
    :param encoding_to_mol: Dictionary mapping molecule encodings to RDKit molecules
    :type encoding_to_mol: Dict[int, Chem.Mol]
    :param encoding_to_rxn: Dictionary mapping reaction encodings to Reaction objects
    :type encoding_to_rxn: Dict[int, Reaction]
    :param reaction_graph: Dictionary representing the reaction graph
    :type reaction_graph: Dict[int, Dict[int, List[int]]]
    :param sequencing_rules: List of sequencing rules
    :type sequencing_rules: List[Motif]
    :param logger: Logger object
    :type logger: Optional[logging.Logger]
    :return: A tuple containing the following:
        - A dictionary mapping molecule encodings to RDKit molecules.
        - A dictionary mapping reaction encodings to Reaction objects.
        - A dictionary representing the reaction graph.
    """
    # retrieve parent molecule
    parent = encoding_to_mol[mol_encoding]
    tag_to_idx = {atom.GetIsotope(): atom.GetIdx() for atom in parent.GetAtoms() if atom.GetIsotope() != 0}

    # create a graph representation of the molecule
    mol_repr = nx.Graph()
    for atom in parent.GetAtoms():
        mol_repr.add_node(atom.GetIdx())
    for bond in parent.GetBonds():
        mol_repr.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())

    # quick filter to see if any sequencing rule can be applied to parent
    if not any([rule(parent) for rule in sequencing_rules]):
        if logger: logger.debug(f"no sequencing rule can be applied to {Chem.MolToSmiles(parent)}")
        return encoding_to_mol, encoding_to_rxn, reaction_graph
    
    # report parent molecule
    if logger: logger.debug(f"parent up for for sequencing: {mol_encoding} {Chem.MolToSmiles(parent)}")

    # a parent given to this mol should only have one chain, that is the assumption, exhaustively apply all chaining rules
    chain = ([], deepcopy(parent))  # current chain, rest
    unfinished_chains = [chain]
    finished_chains = []
    while unfinished_chains:
        current_chain, rest_to_process = unfinished_chains.pop(0)

        new_chains = []
        for rxn in sequencing_rules:
            results = rxn(rest_to_process, logger=logger)
            for result in results:
                if len(result) != 2: 
                    raise ValueError("sequence rule produced more than two products")
                
                # the left product is always the one that is connected to the chain
                # the right product is always the one that is disconnected from the chain
                product1 = result[0]
                product2 = result[1]
                new_chain1 = (current_chain + [product2], product1)

                # get tags of new_chain. Does it form a single connected component
                temp_mol_repr = deepcopy(mol_repr)

                # get all tags from new_chain
                new_chain_tags = set()
                for unit in new_chain1[0]:
                    for atom in unit.GetAtoms():
                        if atom.GetIsotope() != 0:
                            new_chain_tags.add(atom.GetIsotope())

                # convert tags to indices
                new_chain_indices = set([tag_to_idx[tag] for tag in new_chain_tags])

                # remove all atoms that are not in new_chain_indices
                temp_mol_repr.remove_nodes_from([node for node in temp_mol_repr.nodes if node not in new_chain_indices])

                # check if it is a single connected component, this makes sure that the chain is not broken
                if not nx.is_connected(temp_mol_repr):
                    if logger: logger.debug(f"sequence rule {rxn.name} broke the chain")
                    continue

                new_chains.append(new_chain1)

        # if no new chains, add current chain to finished chains
        if not new_chains:
            # add remainder as motif unit (actually starter) to finished chains
            finished_chains.append(current_chain + [rest_to_process])
        else:
            # add new chains to unfinished chains
            unfinished_chains.extend(new_chains)

    # if there are no finished chains, return
    if not finished_chains:
        if logger: logger.debug(f"no chains found for molecule {mol_encoding}")
        return encoding_to_mol, encoding_to_rxn, reaction_graph

    # check what is longest chain and only keep ones that have this length
    max_chain_length = max([len(chain) for chain in finished_chains])
    finished_chains = [chain for chain in finished_chains if len(chain) == max_chain_length]

    # add chaining reaction node to reaction graph, and add children
    for finished_chain in finished_chains:
        num_rxn_nodes = max([enc for enc in encoding_to_rxn]) + 1 if encoding_to_rxn else 1
        reaction_graph[mol_encoding][num_rxn_nodes] = set()

        # get sequence identities and construct sequence representation
        sequence_identities = []
        for sequence_motif in finished_chain:
            identity = match_mol(sequence_motif, logger)
            if identity is not None:
                sequence_identities.append(identity)
            else:
                sequence_identities.append("UNKNOWN")
        sequence_repr = ">".join(sequence_identities)
        if logger: logger.debug(f"sequence: {sequence_repr}")

        # add children to reaction graph
        encoding_to_rxn[num_rxn_nodes] = set([encode_mol(motif) for motif in finished_chain])

    # return modified data structures
    return encoding_to_mol, encoding_to_rxn, reaction_graph


def run_retromol(name: str, smiles: str, logger: Optional[logging.Logger] = None) -> None:
    """Parse a SMILES string.
    
    :param name: Name of the molecule
    :type name: str
    :param smiles: SMILES string
    :type smiles: str
    :param logger: Logger object
    :type logger: Optional[logging.Logger]
    :return: None
    """
    # retrieve reaction rules
    if logger: logger.debug("loading default reaction rules")
    linearization_rules = get_default_linearization_rules()
    if logger: logger.debug(f"loaded {len(linearization_rules)} linearization reaction rules")
    sequencing_rules = get_default_sequencing_rules()
    if logger: logger.debug(f"loaded {len(sequencing_rules)} sequencing reaction rules")

    # preprocess the molecule
    if logger: logger.debug(f"preprocessing molecule {name} with linearization rules")
    encoding_to_mol, encoding_to_rxn, reaction_graph = preprocess_mol(smiles, linearization_rules, logger)
    leaf_nodes = [parent for parent, rxns in reaction_graph.items() if not rxns]
    if logger: logger.debug(f"found {len(leaf_nodes)} leaf nodes")

    # try to sequence all the leaf nodes
    for leaf_node in leaf_nodes:
        mol_encoding = leaf_node
        encoding_to_mol, encoding_to_rxn, reaction_graph = sequence_mol(mol_encoding, encoding_to_mol, encoding_to_rxn, reaction_graph, sequencing_rules, logger)

    return None