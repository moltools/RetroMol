from copy import deepcopy
from typing import Dict, List, Optional, Tuple
import itertools
import logging

from rdkit import Chem
import networkx as nx

from retromol.chem import Reactant, encode_mol, neutralize_mol
from retromol.chem import Reaction


def preprocess_mol(
    reactant: Reactant, 
    reaction_rules: List[Reaction],
    logger: logging.Logger
) -> Tuple[Dict[int, Chem.Mol], Dict[int, Reaction], Dict[int, Dict[int, List[int]]]]:
    """Apply custom rules to linearize a SMILES string.
    
    :param reactant: Reactant object
    :type reactant: Reactant
    :param reaction_rules: List of reaction rules
    :type reaction_rules: List[Reaction]
    :param logger: Logger object
    :type logger: logging.Logger
    :return: A tuple containing the following:
        - A dictionary mapping molecule encodings to RDKit molecules.
        - A dictionary mapping reaction encodings to Reaction objects.
        - A dictionary representing the reaction graph.
    """
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

    # we stored atom indices as isotope numbers, since those persist through reactions
    # retrieve original isotope numbers as tags
    original_taken_tags = reactant.get_tags()

    # create data structures
    encoding_to_mol = {}  # maps mol hash to Chem.Mol
    num_rxn_nodes = 0  # keeps track of number of reaction nodes in reaction graph
    encoding_to_rxn = {}  # maps rxn index to Reaction
    reaction_graph = {}  # as {reactant_mol_encoding: {reaction_encoding: [child_mol_encodings], ...}, ...}

    # create queue
    mols = [deepcopy(reactant.mol)]

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

            # check if reaction_rule contains any ring matching conditions like ;R or ;!R, always make these reactions contested
            if any([reaction_rule.reaction_smarts.find(f";{ring_condition}") != -1 for ring_condition in ["R", "!R"]]):
                contested.append((reaction_rule, match_set))
                continue

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
                    while tag in original_taken_tags or tag in temp_taken_tags:
                        tag += 1
                    atom.SetIsotope(tag)
                    temp_taken_tags.append(tag)
            
            # make sure all atoms have a unique tag
            tagged_atoms = set([atom.GetIsotope() for atom in parent.GetAtoms() if atom.GetIsotope() != 0])
            if len(tagged_atoms) != len(parent.GetAtoms()):
                raise ValueError("not all atoms have a unique tag before applying uncontested reactions")

            # map tags to atomic nums so we can use them to mask atoms
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
                        while tag in temp_taken_tags_uncontested or tag in original_taken_tags or tag in temp_taken_tags:
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


def sequence_mol(
    mol: Chem.Mol,
    sequencing_rules: List[Reaction],
    logger: Optional[logging.Logger] = None
) -> List[List[Chem.Mol]]:
    """Apply custom rules to sequence a molecule into motifs.
    
    :param mol: RDKit molecule
    :type mol: Chem.Mol
    :param sequencing_rules: List of sequencing rules
    :type sequencing_rules: List[Motif]
    :param logger: Logger object
    :type logger: Optional[logging.Logger]
    :return: A list of lists of RDKit molecules, each list representing a motif code
    :rtype: List[List[Chem.Mol]]
    """
    # data structure for saving found motif codes
    finished_chains = []

    # retrieve parent molecule
    tag_to_idx = {atom.GetIsotope(): atom.GetIdx() for atom in mol.GetAtoms() if atom.GetIsotope() != 0}

    # create a graph representation of the molecule
    mol_repr = nx.Graph()
    for atom in mol.GetAtoms():
        mol_repr.add_node(atom.GetIdx())
    for bond in mol.GetBonds():
        mol_repr.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())

    # quick filter to see if any sequencing rule can be applied to parent
    if not any([rule(mol) for rule in sequencing_rules]):
        if logger: logger.debug(f"no sequencing rule can be applied to {Chem.MolToSmiles(mol)}")
        return finished_chains
    
    # report parent molecule
    if logger: logger.debug(f"parent up for for sequencing {Chem.MolToSmiles(mol)}")

    # a parent given to this mol should only have one chain, that is the assumption, exhaustively apply all chaining rules
    chain = ([], deepcopy(mol))  # current chain, rest
    unfinished_chains = [chain]
    while unfinished_chains:
        current_chain, rest_to_process = unfinished_chains.pop(0)

        new_chains = []
        for rxn in sequencing_rules:

            # apply reaction rule
            results = rxn(rest_to_process, logger=logger)

            for result in results:
                if len(result) < 2: 
                    raise ValueError("sequence rule produced less than two products")
                
                if len(result) > 2:
                    if logger: logger.warning(f"sequence rule {rxn.name} produced more than two products, only the first two are considered")
                
                # the left product is always the one that is connected to the chain
                # the right product is always the one that is disconnected from the chain
                product1 = result[0]
                product2 = result[1]
                new_chain1 = (current_chain + [product2], product1)

                # the new found motif has to to match polyketide, alpha amino acid, or beta amino acid backbone smarts
                polyketide_smarts = Chem.MolFromSmarts("[OH]SC~CC(=O)[OH]")
                alpha_amino_acid_smarts = Chem.MolFromSmarts("[NH2]CC(=O)[OH]")
                beta_amino_acid_smarts = Chem.MolFromSmarts("[NH2]CCC(=O)[OH]")

                if (
                    not product2.HasSubstructMatch(polyketide_smarts) and
                    not product2.HasSubstructMatch(alpha_amino_acid_smarts) and
                    not product2.HasSubstructMatch(beta_amino_acid_smarts)
                ):
                    if logger: logger.debug(f"mined motif does not match polyketide, alpha amino acid, or beta amino acid backbone pattern")
                    continue

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
        if logger: logger.debug(f"no chains found for molecule {mol}")
        return finished_chains

    # check what is longest chain and only keep ones that have this length
    max_chain_length = max([len(chain) for chain in finished_chains])
    finished_chains = [chain for chain in finished_chains if len(chain) == max_chain_length]

    if logger:
        logger.debug(f"found {len(finished_chains)} chains for molecule {Chem.MolToSmiles(mol)}")

    return finished_chains
