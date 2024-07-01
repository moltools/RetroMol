# -*- coding: utf-8 -*-

"""Functions for interaction witn PARAS web service."""

import json
import logging
import typing as ty
from enum import Enum

import requests


class ParasInputType(Enum):
    """Input type for PARAS web service.
    
    :cvar FASTA: Fasta format.
    :cvar GENBANK: Genbank format.
    """
    FASTA = "fasta"
    GENBANK = "genbank"


class ParasSubstrateChoice(Enum):
    """Substrate choice for PARAS web service.
    
    :cvar ALL_SUBSTRATES: All substrates.
    :cvar COMMON_SUBSTRATES: Common substrates.
    """
    ALL_SUBSTRATES = "allSubstrates"
    COMMON_SUBSTRATES = "commonSubstrates"


def predict_specificity(
    src: str,
    selected_input_type: ParasInputType = ParasInputType.FASTA,
    save_active_site_signatures: bool = False,
    save_extended_signatures: bool = False,
    save_adenylation_domain_sequences: bool = False,
    selected_substrate_choice: ParasSubstrateChoice = ParasSubstrateChoice.ALL_SUBSTRATES,
    num_predictions_to_report: int = 1,
    use_structure_guided_profile_alignment: bool = False,
    first_separator: str = "|",
    second_separator: str = "_",
    third_separator: str = "-",
) -> ty.List[ty.Tuple[str, ty.List[ty.Tuple[str, float]]]]:
    """Predict specificity of a given sequence.

    :param src: Sequence to predict specificity for.
    :type src: str
    :param selected_input_type: Input type for the sequence, defaults to ParasInputType.FASTA
    :type selected_input_type: ParasInputType, optional
    :param save_active_site_signatures: Save active site signatures, defaults to False
    :type save_active_site_signatures: bool, optional
    :param save_extended_signatures: Save extended signatures, defaults to False
    :type save_extended_signatures: bool, optional
    :param save_adenylation_domain_sequences: Save adenylation domain sequences, defaults to False
    :type save_adenylation_domain_sequences: bool, optional
    :param selected_substrate_choice: Substrate choice, defaults to ParasSubstrateChoice.ALL_SUBSTRATES
    :type selected_substrate_choice: ParasSubstrateChoice, optional
    :param num_predictions_to_report: Number of predictions to report, defaults to 1
    :type num_predictions_to_report: int, optional
    :param use_structure_guided_profile_alignment: Use structure guided profile alignment, defaults to False
    :type use_structure_guided_profile_alignment: bool, optional
    :param first_separator: First separator, defaults to "|"
    :type first_separator: str, optional
    :param second_separator: Second separator, defaults to "_"
    :type second_separator: str, optional
    :param third_separator: Third separator, defaults to "-"
    :type third_separator: str, optional
    :return: predictions
    :rtype: ty.Dict[str, ty.List[ty.Tuple[str, float]]]
    """
    logger = logging.getLogger(__name__)

    url = "https://paras.bioinformatics.nl/api/submit_paras"
    data = {
        "data": {
            "src": src,
            "selectedInputType": selected_input_type.value,
            "saveActiveSiteSignatures": save_active_site_signatures,
            "saveExtendedSignatures": save_extended_signatures,
            "saveAdenylationDomainSequences": save_adenylation_domain_sequences,
            "selectedSubstrateChoice": selected_substrate_choice.value,
            "numPredictionsToReport": num_predictions_to_report,
            "useStructureGuidedProfileAlignment": use_structure_guided_profile_alignment,
            "firstSeparator": first_separator,
            "secondSeparator": second_separator,
            "thirdSeparator": third_separator,
        }
    }
    headers = {"Content-Type": "application/json"}
    response = requests.post(url, data=json.dumps(data), headers=headers).json()
    predictions = []

    if payload := response.get("payload", None):
        if results := payload.get("results", None):
            
            for domain in results:
                domain_id = domain["domain_id"]
                domain_predictions = domain["data"]["predictions"]
                predictions.append((domain_id, domain_predictions))

    return predictions


def label_to_pubchem_cid(label: str) -> int:
    """Convert label to PubChem CID.
    
    :param label: Label to convert.
    :type label: str
    :return: PubChem CID.
    :rtype: int
    """
    return {
        "tryptophan": 6305,
        "arginine": 6322,
        "tyrosine": 6057,
        "glutamine": 5961,
        "leucine": 6106,
        "lysine": 5962,
        "phenylalanine": 6140,
        "R-beta-hydroxytyrosine": "Any",
        "N6-hydroxylysine": 439370,
        "histidine": 6274,
        "glycine": 750,
        "alanine": 5950,
        "N5-hydroxyornithine": 169671,
        "valine": 6287,
        "ornithine": 6262,
        "N5-formyl-N5-hydroxyornithine": 45479225,
        "N5-acetyl-N5-hydroxyornithine": 216641,
        "3-hydroxytyrosine": 836,
        "isoleucine": 6306,
        "beta-hydroxytyrosine": 13309270,
        "4-hydroxyphenylglycine": 92143,
        "3,5-dihydroxyphenylglycine": 108001,
        "2,4-diaminobutyric acid": 470,
        "glutamic acid": 33032,
        "trans-2-crotylglycine": "Any",
        "beta-alanine": 239,
        "serine": 5951,
        "phenylglycine": "Any",
        "methionine": 6137,
        "kynurenine": 846,
        "homotyrosine": 15160483,
        "enduracididine": 15284838,
        "cysteine": 5862,
        "asparagine": 6267,
        "anthranilic acid": 227,
        "O-methyltyrosine": 97118,
        "D-lysergic acid": 6717,
        "6-chlorotryptophan": 65259,
        "4R-propylproline": "Any",
        "3-nitrotyrosine": 235719,
        "2,3-dihydroxybenzoic acid": 19,
        "(S,E)-2-amino-4-decenoic acid": "Any",
        "valinol": 79019,
        "threonine": 6288,
        "succinyl-hydrazinoacetic acid": "Any",
        "succinic semialdehyde": 1112,
        "salicylic acid": 338,
        "quinoxaline-2-carboxylic acid": 96695,
        "pyruvic acid": 1060,
        "pyrrole-2-carboxylic acid": 12473,
        "proline": 145742,
        "piperazic acid": "Any",
        "pipecolic acid": 849,
        "phenylpyruvic acid": 997,
        "phenyllactic acid": 1303,
        "phenazine-1,6-dicarboxylic acid": 193025,
        "para-aminobenzoic acid": 978,
        "p-hydroxymandelate": 328,
        "norcoronamic acid": "Any",
        "nicotinic acid": 938,
        "meta-tyrosine": 6950578,
        "malic acid": 525,
        "linoleic acid": 5280450,
        "lactic acid": 612,
        "homoserine": 779,
        "guanidinoacetic acid": 763,
        "graminine": "Any",
        "fatty acid": "Any",
        "dimethylsulfoniopropionic acid": "Any",
        "dehydrovaline": 15383277,
        "dehydrotryptophan": 129661602,
        "dehydrophenylalanine": 5702627,
        "cysteine branched": "Any",
        "cysteic acid": 25701,
        "coumaric acid": "Any",
        "citrulline": 9750,
        "cinnamic acid": 444539,
        "capreomycidine": "Any",
        "butyric acid": 264,
        "betaine": 247,
        "beta-lysine": 392,
        "beta-hydroxyphenylalanine": 94134,
        "beta-hydroxy-3-hydroxy-O-methyl-5-methyltyrosine": "Any",
        "benzoxazolinate": 5281925,
        "azetidine-2-carboxylic acid": 17288,
        "aspartic acid branched": "Any",
        "aspartic acid": 5960,
        "allo-threonine": 99289,
        "allo-isoleucine": 99288,
        "alaninol": 5126,
        "acetic acid": 176,
        "S-beta-hydroxyenduracididine": "Any",
        "S-beta-hydroxycyclohex-2S-enylalanine": "Any",
        "S-adenosylmethionine": 34755,
        "R-beta-tyrosine": 6934228,  # Manually annotated
        "R-beta-phenylalanine": 6921434,  # Manually annotated
        "R-beta-methyltryptophan": "Any",
        "R-beta-methylphenylalanine": "Any",
        "R-beta-hydroxyphenylalanine": "Any",
        "R-aza-beta-tyrosine": "Any",
        "R-3-hydroxy-3-methylproline": "Any",
        "N5-trans-anhydromevalonyl-N5-hydroxyornithine": "Any",
        "N5-cis-anhydromevalonyl-N5-hydroxyornithine": "Any",
        "N5-acetyl-hydroxyornithine": "Any",
        "N1-methoxytryptophan": "Any",
        "N-hydroxyvaline": 22859864,
        "N-formylglycine": 75606,
        "D-phenyllactic acid": 444718,
        "D-leucine": 439524,
        "D-glutamic acid branched": "Any",
        "D-aspartic acid branched": "Any",
        "D-alanine": 71080,
        "Compound 4 (formed by the decarboxylative condensation of L-Phe and succinyl-CoA)": "Any",
        "An acid hydrazine polyene (intermediate 14)": "Any",
        "6S-methyl-pipecolic acid": "Any",
        "6-methylsalicylic acid": 11279,
        "6-hydroxy-tetrahydro-isoquinoline-3-carboxylic acid": "Any",
        "6-chloro-4-hydroxyindole-3-carboxylic acid": "Any",
        "6-chloro-4-hydroxy-1-methyl-indole-3-carboxylic acid": "Any",
        "6,7-dichlorotryptophan": 91530202,
        "5S-methylproline": "Any",
        "5-methylorsellinic acid": 71365885,
        "5-methoxytyrosine": "Any",
        "5-chlorotryptophan": 644330,
        "5-chloroanthranilic acid": 12476,
        "5-aminolevulinic acid": 137,
        "5,5-dimethylpipecolic acid": "Any",
        "4S-propenylproline": "Any",
        "4S-methylproline": "Any",
        "4S-methylazetidine-2S-carboxylic acid": "Any",
        "4S-hydroxylysine": "Any",
        "4S-acetyl-5S-methylproline": "Any",
        "4S,5-dihydroxy-2S-aminopentanoic acid": "Any",
        "4R-methylproline": 57377078,
        "4R-hydroxyproline": "Any",
        "4R-E-butenyl-4R-methylthreonine": "Any",
        "4-nitrotryptophan": 90403874,
        "4-methylproline": 352051,
        "4-methoxytryptophan": "Any",
        "4-hydroxyphenylpyruvic acid": 979,
        "4-hydroxyglutamine": 54409926,
        "4-hydroxybenzoic acid": 135,
        "4-hydroxy-D-kynurenine": "Any",
        "4-hydroxy-3-nitrobenzoic acid": 12033,
        "4-aminobutyric acid": 119,
        "4-amino-2-hydroxy-3-isopropoxybenzoic acid": 118989950,
        "4-acetamidopyrrole-2-carboxylic acid": "Any",
        "4,5-dehydroarginine": "Any",
        "3S-methylproline": "Any",
        "3S-methylleucine": "Any",
        "3S-methylaspartic acid branched": "Any",
        "3S-methylaspartic acid": "Any",
        "3S-methyl-D-aspartic acid branched": "Any",
        "3S-hydroxyleucine": "Any",
        "3S-hydroxyasparagine": "Any",
        "3S-hydroxy-4S-methylproline": "Any",
        "3S-cyclohex-2-enylalanine": "Any",
        "3R-methylglutamic acid": 439377, # Manually annotated
        "3R-methylbeta-alanine": "Any",
        "3R-hydroxyhomotyrosine": "Any",
        "3R-hydroxyaspartic acid": "Any",
        "3R-hydroxyasparagine": "Any",
        "3R-hydroxy-2,4-diaminobutyric acid": "Any",
        "3-methylglutamic acid": 237657,
        "3-methylasparagine": "Any",
        "3-methoxyaspartic acid": "Any",
        "3-methoxyanthranilic acid": 255720,
        "3-hydroxyvaline": 192763,
        "3-hydroxyquinaldic acid": 194859,
        "3-hydroxypicolinic acid": 13401,
        "3-hydroxyleucine": 277776,
        "3-hydroxykynurenine": 89,
        "3-hydroxyglutamine": 22855697,
        "3-hydroxyaspartic acid": 5425,
        "3-hydroxyasparagine": 152191,
        "3-hydroxy-para-aminobenzoic acid": "Any",
        "3-hydroxy-O-methyl-5-methyltyrosine": "Any",
        "3-hydroxy-4-methylproline": 53731613,
        "3-aminoisobutyric acid": 64956,
        "3-amino-6-hydroxy-2-piperidone": 54082380,
        "3-amino-4-hydroxybenzoic acid": 65083,
        "3-amino-2,4-dihydroxybenzoic acid": 24899360,
        "3-(3-pyridyl)-alanine": "Any",
        "3-(2-nitrocyclopropylalanine)": "Any",
        "3,4-dihydroxybenzoic acid": 72,
        "3,4-dehydrolysine": 129689982,
        "2S-methyl-3-oxobutyrine": "Any",
        "2S-hydroxyisovaleric acid": "Any",
        "2S-hydroxyisocaproic acid": "Any",
        "2S-aminooctanoic acid": "Any",
        "2S-aminododecanoic acid": "Any",
        "2S-aminodecanoic acid": "Any",
        "2S-amino-9,10-epoxy-8-oxodecanoic acid": "Any",
        "2S-amino-8-oxodecanoic acid": "Any",
        "2S,3S-diaminobutyric acid": "Any",
        "2R-hydroxyisovaleric acid": "Any",
        "2R-hydroxy-3-methylpentanoic acid": "Any",
        "2-sulfamoylacetic acid": 11083961,
        "2-methylserine": 94309,
        "2-ketoisovaleric acid": 49,
        "2-ketoisocaproic acid": 70,
        "2-ketoglutaric acid": 51,
        "2-chlorobenzoic acid": 8374,
        "2-carboxy-6-hydroxyoctahydroindole": "Any",
        "2-aminoisobutyric acid": 6119,
        "2-aminobutyric acid": 6657,
        "2-aminoadipic acid": 469,
        "2-amino-6-hydroxy-4-methyl-8-oxodecanoic acid": "Any",
        "2-amino-3-hydroxycyclopent-2-enone": 11073376,
        "2-amino-3,5-dimethyl-4-hexenoic Acid": 19084654,
        "2-(1-methylcyclopropyl)-D-glycine": "Any",
        "2,4-dihydroxypentanoic acid": 13644067,
        "2,3-dihydroxyhexadecanoic acid": 19591836,
        "2,3-dihydroxy-para-aminobenzoic acid": "Any",
        "2,3-diaminopropionic acid": 364,
        "2,3-diaminobutyric acid": 5316628,
        "1-pyrroline-5-carboxylic acid": 1196,
        "1-aminocyclopropane-1-carboxylic acid": 535,
        "1-(1,1-dimethylallyl)-tryptophan": "Any",
        "(E)-4-methylhex-2-enoic acid": 5463168,
        "(4S)-5,5,5-trichloroleucine": "Any",
        "(2S,6R)-diamino-(5R,7)-dihydroxy-heptanoic acid": 121514024,
        "(2S,3R)-2-amino-3-hydroxy-4-(4-nitrophenyl)butanoic acid": "Any",
    }.get(label, None)
