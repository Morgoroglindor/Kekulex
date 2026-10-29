"""
templates.py - ready-made structures for the Templates menu, as SMILES (laid out by RDKit when inserted).

    TEMPLATES = {category: [(name, smiles)]}
"""
TEMPLATES = {
    'Rings and ring systems': [
        ('Cyclopropane', 'C1CC1'), ('Cyclobutane', 'C1CCC1'), ('Cyclopentane', 'C1CCCC1'), ('Cyclohexane', 'C1CCCCC1'),
        ('Cycloheptane', 'C1CCCCCC1'), ('Benzene', 'c1ccccc1'), ('Naphthalene', 'c1ccc2ccccc2c1'),
        ('Anthracene', 'c1ccc2cc3ccccc3cc2c1'), ('Phenanthrene', 'c1ccc2c(c1)ccc1ccccc12'),
        ('Biphenyl', 'c1ccc(-c2ccccc2)cc1'), ('Indane', 'c1ccc2c(c1)CCC2'), ('Tetralin', 'c1ccc2c(c1)CCCC2'),
        ('Decalin', 'C1CCC2CCCCC2C1'), ('Norbornane', 'C1CC2CCC1C2'), ('Adamantane', 'C1C2CC3CC1CC(C2)C3'),
        ('Spiro[4.5]decane', 'C1CCC2(C1)CCCCC2'), ('Steroid skeleton', 'C1CCC2C(C1)CCC1C2CCC2CCCC21'),
    ],
    'Heterocycles': [
        ('Pyridine', 'c1ccncc1'), ('Pyrimidine', 'c1cncnc1'), ('Pyrazine', 'c1cnccn1'), ('Pyridazine', 'c1ccnnc1'),
        ('1,3,5-Triazine', 'c1ncncn1'), ('Pyrrole', 'c1cc[nH]c1'), ('Furan', 'c1ccoc1'), ('Thiophene', 'c1ccsc1'),
        ('Imidazole', 'c1cnc[nH]1'), ('Pyrazole', 'c1cn[nH]c1'), ('Oxazole', 'c1cocn1'), ('Isoxazole', 'c1cnoc1'),
        ('Thiazole', 'c1cscn1'), ('1,2,3-Triazole', 'c1c[nH]nn1'), ('Tetrazole', 'c1nnn[nH]1'),
        ('Indole', 'c1ccc2[nH]ccc2c1'), ('Benzofuran', 'c1ccc2occc2c1'), ('Benzothiophene', 'c1ccc2sccc2c1'),
        ('Quinoline', 'c1ccc2ncccc2c1'), ('Isoquinoline', 'c1ccc2cnccc2c1'), ('Quinoxaline', 'c1ccc2nccnc2c1'),
        ('Benzimidazole', 'c1ccc2[nH]cnc2c1'), ('Benzoxazole', 'c1ccc2ocnc2c1'), ('Purine', 'c1ncc2[nH]cnc2n1'),
        ('Carbazole', 'c1ccc2c(c1)[nH]c1ccccc12'), ('Acridine', 'c1ccc2nc3ccccc3cc2c1'),
        ('Piperidine', 'C1CCNCC1'), ('Piperazine', 'C1CNCCN1'), ('Morpholine', 'C1COCCN1'), ('Pyrrolidine', 'C1CCNC1'),
        ('Azetidine', 'C1CNC1'), ('Aziridine', 'C1CN1'), ('Oxirane', 'C1CO1'), ('Tetrahydrofuran', 'C1CCOC1'),
        ('Tetrahydropyran', 'C1CCOCC1'), ('1,4-Dioxane', 'C1COCCO1'), ('Coumarin', 'O=c1ccc2ccccc2o1'),
        ('Chromone', 'O=c1ccoc2ccccc12'), ('Flavone', 'O=c1cc(-c2ccccc2)oc2ccccc12'),
        ('Porphine', 'c1cc2cc3ccc(cc4ccc(cc5ccc(cc1n2)[nH]5)n4)[nH]3'),
    ],
    'Functional groups': [
        ('Carboxylic acid', 'CC(=O)O'), ('Ester', 'CC(=O)OC'), ('Amide', 'CC(N)=O'), ('Acid chloride', 'CC(=O)Cl'),
        ('Anhydride', 'CC(=O)OC(C)=O'), ('Aldehyde', 'CC=O'), ('Ketone', 'CC(C)=O'), ('Nitrile', 'CC#N'),
        ('Imine', 'CC=NC'), ('Oxime', 'CC=NO'), ('Hydrazone', 'CC=NN'), ('Enamine', 'C=CN(C)C'), ('Enol ether', 'C=COC'),
        ('Nitro', 'C[N+](=O)[O-]'), ('Azide', 'CN=[N+]=[N-]'), ('Diazo', 'C=[N+]=[N-]'), ('Isocyanate', 'CN=C=O'),
        ('Isothiocyanate', 'CN=C=S'), ('Carbamate', 'CNC(=O)OC'), ('Urea', 'CNC(=O)NC'), ('Carbonate', 'COC(=O)OC'),
        ('Thiol', 'CS'), ('Thioether', 'CSC'), ('Disulfide', 'CSSC'), ('Sulfoxide', 'CS(C)=O'), ('Sulfone', 'CS(C)(=O)=O'),
        ('Sulfonamide', 'CS(=O)(=O)N'), ('Sulfonate ester', 'CS(=O)(=O)OC'), ('Phosphate ester', 'COP(=O)(OC)OC'),
        ('Boronic acid', 'CB(O)O'), ('Epoxide', 'CC1CO1'), ('Acetal', 'CC(OC)OC'), ('Allene', 'C=C=C'),
    ],
    'Reagents': [
        ('NBS', 'O=C1CCC(=O)N1Br'), ('NCS', 'O=C1CCC(=O)N1Cl'), ('DMAP', 'CN(C)c1ccncc1'), ('DBU', 'C1CCC2=NCCCN2CC1'),
        ('DIPEA', 'CCN(C(C)C)C(C)C'), ('Triethylamine', 'CCN(CC)CC'), ('Pyridine', 'c1ccncc1'), ('Imidazole', 'c1cnc[nH]1'),
        ('DCC', 'C(=NC1CCCCC1)=NC1CCCCC1'), ('EDC', 'CCN=C=NCCCN(C)C'), ('HOBt', 'On1nnc2ccccc21'),
        ('CDI', 'O=C(n1ccnc1)n1ccnc1'), ('Boc2O', 'CC(C)(C)OC(=O)OC(=O)OC(C)(C)C'),
        ('TsCl', 'Cc1ccc(S(=O)(=O)Cl)cc1'), ('MsCl', 'CS(=O)(=O)Cl'), ('Tf2O', 'O=S(=O)(OS(=O)(=O)C(F)(F)F)C(F)(F)F'),
        ('TBSCl', 'CC(C)(C)[Si](C)(C)Cl'), ('TMSCl', 'C[Si](C)(C)Cl'), ('Acetic anhydride', 'CC(=O)OC(C)=O'),
        ('Thionyl chloride', 'O=S(Cl)Cl'), ('Oxalyl chloride', 'O=C(Cl)C(=O)Cl'), ('POCl3', 'O=P(Cl)(Cl)Cl'),
        ('mCPBA', 'OOC(=O)c1cccc(Cl)c1'), ('TEMPO', 'CC1(C)CCCC(C)(C)N1[O]'), ('DMP (Dess–Martin)', 'CC(=O)OI1(OC(C)=O)(OC(C)=O)OC(=O)c2ccccc21'),
        ('NaBH4', '[Na+].[BH4-]'), ('LiAlH4', '[Li+].[AlH4-]'), ('DIBAL-H', 'CC(C)C[AlH]CC(C)C'),
        ('LDA', '[Li+].CC(C)[N-]C(C)C'), ('n-BuLi', '[Li]CCCC'), ('NaH', '[Na+].[H-]'), ('KOtBu', '[K+].CC(C)(C)[O-]'),
        ('Triphenylphosphine', 'c1ccc(P(c2ccccc2)c2ccccc2)cc1'), ('DEAD', 'CCOC(=O)N=NC(=O)OCC'),
        ('Lawesson\'s reagent', 'COc1ccc(P2(=S)SP(=S)(c3ccc(OC)cc3)S2)cc1'),
    ],
    'Solvents': [
        ('Water', 'O'), ('Methanol', 'CO'), ('Ethanol', 'CCO'), ('Isopropanol', 'CC(C)O'), ('Acetonitrile', 'CC#N'),
        ('Acetone', 'CC(C)=O'), ('DMSO', 'CS(C)=O'), ('DMF', 'CN(C)C=O'), ('Dichloromethane', 'ClCCl'),
        ('Chloroform', 'ClC(Cl)Cl'), ('Ethyl acetate', 'CCOC(C)=O'), ('Diethyl ether', 'CCOCC'), ('THF', 'C1CCOC1'),
        ('1,4-Dioxane', 'C1COCCO1'), ('Toluene', 'Cc1ccccc1'), ('Hexane', 'CCCCCC'), ('Formic acid', 'O=CO'),
        ('Acetic acid', 'CC(=O)O'), ('Trifluoroacetic acid', 'OC(=O)C(F)(F)F'), ('Sulfuric acid', 'O=S(=O)(O)O'),
        ('Nitric acid', 'O=[N+]([O-])O'),
    ],
}
