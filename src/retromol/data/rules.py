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
        "reaction_smarts": r"[N:1][CH3:2]>>[N:1].[C:2]"
    },
    {
        "name": "s-methylation",
        "reaction_smarts": r"[S:1][CH3:2]>>[S:1].[CH4:2]",
    },
    {
        "name": "glycosyltransferase", 
        "reaction_smarts": r"[C:1][O:2][C:3]1[O:4][C:5][C:6][C:7][C:8]1>>[C:1][OH:2].[OH][C:3]1[O:4][C:5][C:6][C:7][C:8]1"
    },
    {
        "name": "etherification",
        "reaction_smarts": r"[*:1]-[CH0:2]1(-[OH:3])-[O:4]-[CH1:5](-[*:6])-[C:7][C:8][C:9]1>>([*:1]-[CH0:2]1=[OH0:3].[OH1:4]-[CH1:5](-[*:6])-[C:7][C:8][C:9]1)",
    },
    {
        "name": "etherification",
        "reaction_smarts": r"[*&D2,D3:3]-[CH1:1]-[C&D3:2]1-[O:4]-[CH1&D3:5](-[*:6])-[C:7][C:8][C:9]1>>([*&D2,D3:3]-[CH0:1]=[C:2]1.[OH:4]-[CH1:5](-[*:6])-[C:7][C:8][C:9]1)",
    },
    {
        "name": "etherification",
        "reaction_smarts": r"[*&D2,D3:3]-[CH2:1]-[C&D3:2]1-[O:4]-[CH1&D3:5](-[*:6])-[C:7][C:8][C:9]1>>([*&D2,D3:3]-[CH1:1]=[C:2]1.[O:4]=[CH0:5](-[*:6])-[C:7][C:8][C:9]1)",
    },
    {
        "name": "etherification",
        "reaction_smarts": r"[*:1]-[CH1:2]([C:3]~[C:4]~[C:5]1)-[O:6]-[CH1:7]1-[CH2:8]-[*:9]>>[*:1]-[C:2](-[OH:6])[C:3]~[C:4]~[C:5][CH0:7]=[CH1:8]-[*:9]",
    },
    {
        "name": "etherification",
        "reaction_smarts": r"[*:1]-[CH1:2]([C:3]~[C:4]~[C:5]1)-[O:6]-[CH1:7]1-[CH1:8]-[*:9]>>[*:1]-[C:2](-[OH:6])[C:3]~[C:4]~[C:5][CH0:7]=[CH0:8]-[*:9]",
    },
    {
        "name": "etherification",
        "reaction_smarts": r"[CH2:2]([C:3]~[C:4]~[C:5]1)-[O:6]-[CH1:7]1-[CH1:8]-[*:9]>>[C:2](-[OH:6])[C:3]~[C:4]~[C:5][CH0:7]=[CH0:8]-[*:9]",
    },
    {
        "name": "etherification",
        "reaction_smarts": r"[CH2:2]([C:3]~[C:4]~[C:5]1)-[O:6]-[CH1:7]1-[CH1:8]-[*:9]>>[C:2](-[OH:6])[C:3]~[C:4]~[C:5][CH0:7]=[CH0:8]-[*:9]",
    },
    {
        "name": "etherification",
        "reaction_smarts": r"[CH2:2]([C:3]~[C:4]~[C:5]1)-[O:6]-[CH1:7]1-[CH2:8]-[*:9]>>[C:2](-[OH:6])[C:3]~[C:4]~[C:5][CH0:7]=[CH1:8]-[*:9]",
    },
    {
        "name": "etherification",
        "reaction_smarts": r"[*&D2,D3:3]-[CH1:1]-[C&D3:2]1-[O:4]-[CH1&D3:5](-[*:6])-[C:7][C:9]1>>([*&D2,D3:3]-[CH0:1]=[C:2]1.[OH:4]-[CH1:5](-[*:6])-[C:7][C:9]1)",
    },
    {
        "name": "macrolactonization", 
        "reaction_smarts": r"[C;R:1][C;R:2](=[O:3])[O;R:4][C;R:5]>>([C:1][C:2](=[O:3])[OH].[OH:4][C:5])"
    },
    {
        "name": "ester bond", 
        "reaction_smarts": r"[C:1][C;!R:2](=[O:3])[O;!R:4][C:5]>>[C:1][C:2](=[O:3])[OH].[OH:4][C:5]"
    },
    {
        "name": "macrolactonethionization", 
        "reaction_smarts": r"[C;R:1][C;R:2](=[O:3])[S;R:4][C;R:5]>>([C:1][C:2](=[O:3])[OH].[SH:4][C:5])"
    },
    {
        "name": "thio-ester bond", 
        "reaction_smarts": r"[C:1][C;!R:2](=[O:3])[S;!R:4][C:5]>>[C:1][C:2](=[O:3])[OH].[SH:4][C:5]"
    },
    # {
    #     "name": "carbocyclization", 
    #     "reaction_smarts": r"[C:1][C:2]1[C:3]([C:4])[C:5]=[C:6][C:7]2[C:8]~[C:9]~[C:10]~[C:11][C:12]21>>([C:1][C:2]1.[C:3]([C:4])[C:5]=[C:6][C:7][C:8]~[C:9]~[C:10]~[C:11]~[C:12]1)"
    # },
    {
        "name": "carbocyclization",
        "reaction_smarts": r"[*:1][C:2]1[C:3]([*:4])[C:5]=[C:6][CH1:7]([CH2;R:8])[C:9]([C;R:10])1>>([*:1][C:2]=1.[C:3]([*:4])[C:5]=[C:6][CH1:7](=[CH1;R:8]).[C:9]([C;R:10])=1)",
    },
    {
        "name": "carbocyclization",
        "reaction_smarts": r"[*:1][C:2]1[C:3]([*:4])[C:5]=[C:6][CH1:7]([CH1;R:8])[C:9]([C;R:10])1>>([*:1][C:2]=1.[C:3]([*:4])[C:5]=[C:6][CH1:7](=[CH0;R:8]).[C:9]([C;R:10])=1)",
    },
    {
        "name": "dihalogenation", 
        "reaction_smarts": r"[C&D4:2]([Cl,Br,I:3])([Cl,Br,I:4])>>[C:2].[Cl,Br,I:3].[Cl,Br,I:4]"
    },
    {
        "name": "carbamic acid", 
        "reaction_smarts": r"[NH2:1]-[CH0:2](=[O:3])-[OH0:4]-[*:5]>>[NH2:1]-[CH0:2](=[O:3])[OH].[OH1:4]-[*:5]"
    },
    {
        "name": "acetic acid", 
        "reaction_smarts": r"[CH3:1]-[CH0:2](=[O:3])-[OH0:4]-[*:5]>>[CH3:1]-[CH0:2](=[O:3])[OH].[OH1:4]-[*:5]"
    },
    {
        "name": "cyanide",
        "reaction_smarts": r"[C:1][C&D2:2]#[N&D1:3]>>[C:1].[C:2]#[N:3]"
    },
    {
        "name": "disulfide bridge",
        "reaction_smarts": r"[C;R:1][C:2][S:3][S:4][C:5][C;R:6]>>([C;R:1][C:2][SH:3].[SH:4][C:5][C;R:6])",
    },
    {
        "name": "disulfide bridge",
        "reaction_smarts": r"[C:1][C:2][S:3][S:4][C:5][C:6]>>[C;R:1][C:2][SH:3].[SH:4][C:5][C;R:6]",
    },
    {
        "name": "hydroxyl sulfonation",
        "reaction_smarts": r"[*:1][S:2](=[O:3])(=[O:4])[OH:5]>>[*:1].[OH][S:2](=[O:3])(=[O:4])[OH:5]",
    },
    {
        "name": "reduction 1",
        "reaction_smarts": r"[NH2:1][CH0:2]([OH:3])[C:4](=[O:5])[OH:6]>>[NH2:1][CH1:2][C:4](=[O:5])[OH:6].[OH2:3]",
    },
    {
        "name": "reduction 2",
        "reaction_smarts": r"[NH2:1][CH0:2]([OH:3])[CH2:4][C:5](=[O:6])[OH:7]>>[NH2:1][CH1:2][CH2:4][C:5](=[O:6])[OH:7].[OH2:3]",
    },
    # {
    #     "name": "reduction 3",
    #     "reaction_smarts": r"[NH1:1]-[CH1:2]-[CH2:3]-[OH:4]>>[NH1:1]-[CH1:2]-[CH0:3](=[O])-[OH:4]",
    # },
    {
        "name": "oxazole",
        "reaction_smarts" : r"[C:1][c:2]1[o:3][c:4][c:5]([C:6])[n:7]1>>([C:1]-[CH0:2]-1(=O).[OH1:3]-[CH2:4]-[CH1:5](-[C:6])-[NH1:7]-1)",
    },
    {
        "name": "oxazoline",
        "reaction_smarts" : r"[C:1][C:2]=1[OH0:3][C:4][C:5]([C:6])[N:7]1>>([C:1]-[CH0:2]-1(=O).[OH1:3]-[CH2:4]-[CH1:5](-[C:6])-[NH1:7]-1)",
    },
    {
        "name": "thiazole",
        "reaction_smarts": r"[C:1][c:2]1[s:3][c:4][c:5][n:7]1>>([C:1]-[CH0:2]-1(=O).[SH1:3]-[CH2:4]-[CH1:5]-[NH1:7]-1)",
    },
    {
        "name": "thiazoline tautomerization",
        "reaction_smarts" : r"[O:1]=[CH0:2]-1[SH0:3][C:4][C:5]([C:6])[NH1:7]1>>[OH:1]-[CH0:2]=1[SH0:3][C:4][C:5]([C:6])[NH0:7]1",
    },
    {
        "name": "thiazoline",
        "reaction_smarts" : r"[C,O:1][C:2]=1[SH0:3][C:4][C:5][N:7]1>>([C,O:1]-[CH0:2]-1(=O).[SH1:3]-[CH2:4]-[C:5]-[NH1:7]-1)",
    },
    {
        "name": "thiazoline",
        "reaction_smarts" : r"[*:1][C:2]([C:3][SH0:4]1)[NH1:5][C:6]1[*:7]>>[*:1][C:2]([C:3][SH1:4])[NH1:5][CH0:6](=[O])[*:7]",
    },
    {
        "name": "tetramate",
        "reaction_smarts" : r"[*:1]-[C:2](-[OH:3])=[C:4]1[C:5](=[O:6])[N:7][C:8][C:9]1=[O:10]>>[*:1]-[C:2](=[OH0:3])-[CH2:4][C:5](=[O:6])[N:7][C:8][C:9](-[OH])=[O:10]",
    },
    {
        "name": "tetramate",
        "reaction_smarts" : r"[*:1]-[C:2](-[OH:3])-[C:4]=1[C:5](=[O:6])[N:7][C:8][C:9]1-[OH:10]>>[*:1]-[C:2](=[OH0:3])-[CH2:4][C:5](=[O:6])[N:7][C:8][C:9](-[OH])=[O:10]",
    },
    {
        "name": "tetronate",
        "reaction_smarts" : r"[*:1]-[C:2](-[OH:3])=[C:4]1[C:5](=[O:6])[O:7][C:8][C:9]1=[O:10]>>[*:1]-[C:2](=[OH0:3])-[CH2:4][C:5](=[O:6])[O:7][C:8][C:9](-[OH])=[O:10]",
    },
    {
        "name": "tetronate",
        "reaction_smarts" : r"[*:1]-[C:2](-[OH:3])-[C:4]=1[C:5](=[O:6])[O:7][C:8][C:9]1-[OH:10]>>[*:1]-[C:2](=[OH0:3])-[CH2:4][C:5](=[O:6])[O:7][C:8][C:9](-[OH])=[O:10]",
    },
    {
        "name": "kirromycin-like",
        "reaction_smarts": r"[c:1]1([C:2](=[O:8])[*:9])[c:3](~[O:10])[c:4][c:5][n:6][c:7](~[O:11])1>>[CH0:3](=[O:10])(-[OH])[CH2:4][CH2:5][NH1:6][CH0:7](=[O:11])[CH2:1]([CH0:2](=[O:8])[*:9])",
    },
    {
        "name": "dioxane etherifcation",
        "reaction_smarts": r"[*:6][C:1]1[C:2][C:3]([*:7])[O:4][CH2][O:5]1>>[*:6][C:1]([OH:5])[C:2][C:3]([OH:4])[*:7]",
    },
    {
        "name": "amido transferase",
        "reaction_smarts": r"[NH1:1][CH1:2]([C:3](=[O:4])[NH2:5])[CH2:6][NH2:7]>>[NH1:1][CH1:2]([C:3](=[O:4])[OH:5])[CH2:6][NH2:7]",
    },
    {
        "name": "hydrogenation",
        "reaction_smarts": r"[CH2:1][CH1:2]([OH1:3])[CH2:4][OH1:5]>>[CH2:1][CH0:2](=[OH0:3])[OH1].[CH3:4][OH1:5]"
    },
    {
        "name": "beta-lactam",
        "reaction_smarts": r"[C:1][N:2][C:3]1[C:4]2[S:5][C:6][C:7]([C:8])[N:9]2[C:10]1=[O:11]>>([C:1][N:2][C:3]1[C:4][S:5].[C:6][C:7]([C:8])[N:9][C:10]1=[O:11])",
    },
    {
        "name": "DAOC-synthetase",
        "reaction_smarts": r"[C:1][N:2][C:3]1[C:4]2[S:5][C:6][CH0:7]=[CH0:8]([C:9])[N:10]2[C:11]1=[O:12]>>[C:1][N:2][C:3]1[C:4]2[S:5][CH0:7]([CH3:6])[CH1:8]([C:9])[N:10]2[C:11]1=[O:12]",
    }
]

SEQUENCING_RULES = [
    {
        "name": "pks (saturated)", 
        "reaction_smarts": r"[C,c:1][C;!R:2]-[C;!R:3]-[C:4](=[O:5])[OH:6]>>[C:1]C(=O)[OH].[OH][S][C:2]-[C:3]-[C:4](=[O:5])[OH:6]"
    },
    {
        "name": "pks (unsaturated)", 
        "reaction_smarts": r"[C,c:1][C;!R:2]=[C;!R:3]-[C:4](=[O:5])[OH:6]>>[C:1]C(=O)[OH].[OH][S][C:2]=[C:3]-[C:4](=[O:5])[OH:6]"
    },
    {
        "name": "pks (unsaturated and shifted)", 
        "reaction_smarts": r"[C,c:1]=[C;!R:2]-[C&D2,C&D3;!R:3]-[C:4](=[O:5])[OH:6]>>[C:1]C(=O)[OH].[OH][S][C:2]=[C:3]-[C:4](=[O:5])[OH:6]"
    },
    {
        "name": "pks (shifted and late stage oxidation)", 
        "reaction_smarts": r"[C,c:1]=[C;!R:2]-[C;!R:3](-[OH:7])-[C:4](=[O:5])[OH:6]>>[C:1]C(=O)[OH].[OH][S][C:2]=[C:3]-[C:4](=[O:5])[OH:6].[O:7]"
    },
    {
        "name": "adenylation domain", 
        "reaction_smarts": r"[*:1][C:2](=[O:3])[NH1:4][C:5][C:6](=[O:7])[OH:8]>>[C:1][C:2](=[O:3])[OH].[NH2:4][C:5][C:6](=[O:7])[OH:8]"
    },
]