# run deorphanization: match compounds to antismash database records -- best matches also match on producing organism (preselect!)

# note that we have many compounds from unknown organisms:
#  organism_type | organism_genus | organism_species | compound_count 
# ---------------+----------------+------------------+----------------
#  bacterium     | Streptomyces   | sp.              |            653
#  fungus        | Unknown-fungus | sp.              |            525
#  fungus        | Penicillium    | sp.              |            216
#  fungus        | Aspergillus    | sp.              |            214
#  fungus        | Ganoderma      | lucidum          |            208
#  bacterium     | Lyngbya        | majuscula        |            171
#  fungus        | Aspergillus    | versicolor       |            159
#  bacterium     | Microcystis    | sp.              |            148
#  fungus        | Aspergillus    | terreus          |            144
#  fungus        | Unknown        | sp.              |            138