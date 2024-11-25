LINEARIZATION_RULES = [
    { 
        "name": "epoxidation", 
        "reaction_smarts": r"[C:1]1[O:2][C:3]1>>[C:1]=[C:3].[O:2]"
    },
    {
        "name": "o-methylation", 
        "reaction_smarts": r"[O&D2:1][CH3:2]>>[O:1].[C:2]"
    },
    {
        "name": "n-methylation",
        "reaction_smarts": r"[N&D2:1][CH3:2]>>[N:1].[C:2]"
    },
    {
        "name": "n-dimethylation", 
        "reaction_smarts": r"[N&D3:1]([CH3:2])[CH3:3]>>[N:1].[C:2].[C:3]"
    },
    {
        "name": "glycosyltransferase", 
        "reaction_smarts": r"[C:1][O:2][C:3]1[O:4][C:5][C:6][C:7][C:8]1>>[C:1][OH:2].[OH][C:3]1[O:4][C:5][C:6][C:7][C:8]1"
    },
    {
        "name": "etherification", 
        "reaction_smarts": r"[C:8][C:1]1[C:2][C:3][C:4][C:5]([O:6]1)(-[OH:7])[C:9]>>[C:8][C:1](-[OH:6])[C:2][C:3][C:4][C:5](=[O:7])[C:9]"
    },
    {
        "name": "macrolactonization", 
        "reaction_smarts": r"[C;R:1][C;R:2](=[O:3])[O;R:4][C;R:5]>>([C:1][C:2](=[O:3])[OH].[OH:4][C:5])"
    },
    {
        "name": "carbocyclization", 
        "reaction_smarts": r"[C:1][C:2]1[C:3]([C:4])[C:5]=[C:6][C:7]2[C:8]~[C:9]~[C:10]~[C:11][C:12]21>>([C:1][C:2]1.[C:3]([C:4])[C:5]=[C:6][C:7][C:8]~[C:9]~[C:10]~[C:11]~[C:12]1)"
    },
    {
        "name": "dihalogenation", 
        "reaction_smarts": r"[C&D4:2]([Cl,Br,I:3])([Cl,Br,I:4])>>[C:2].[Cl,Br,I:3].[Cl,Br,I:4]"
    },
    {
        "name": "carbamic_acid", 
        "reaction_smarts": r"[NH2:1]-[CH0:2](=[O:3])-[OH0:4]-[*:5]>>[NH2:1]-[CH0:2](=[O:3])[OH].[OH1:4]-[*:5]"
    },
]

SEQUENCING_RULES = [
    {
        "name": "pks_saturated", 
        "reaction_smarts": r"[C,c;!R:1][C;!R:2]-[C;!R:3]-[C:4](=[O:5])[OH:6]>>[C:1]C(=O)[OH].[OH][S][C:2]-[C:3]-[C:4](=[O:5])[OH:6]"
    },
    {
        "name": "pks_unsaturated", 
        "reaction_smarts": r"[C,c;!R:1][C;!R:2]=[C;!R:3]-[C:4](=[O:5])[OH:6]>>[C:1]C(=O)[OH].[OH][S][C:2]=[C:3]-[C:4](=[O:5])[OH:6]"
    },
    {
        "name": "pks_unsaturated_shifted", 
        "reaction_smarts": r"[C,c;!R:1]=[C;!R:2]-[C&D2,C&D3;!R:3]-[C:4](=[O:5])[OH:6]>>[C:1]C(=O)[OH].[OH][S][C:2]=[C:3]-[C:4](=[O:5])[OH:6]"
    },
    {
        "name": "adenylation_domain", 
        "reaction_smarts": r"[*:1][C:2](=[O:3])[NH1:4][C:5][C:6](=[O:7])[OH:8]>>[C:1][C:2](=[O:3])[OH].[NH2:4][C:5][C:6](=[O:7])[OH:8]"
    },
]