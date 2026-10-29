"""
templates.py - ready-made structures for the Templates menu, as SMILES (laid out by RDKit when inserted).

    TEMPLATES = {category: [(name, smiles)]}
"""
TEMPLATES = {
    'Rings and ring systems': [
        ('Cyclopropane', 'C1CC1'), ('Cyclobutane', 'C1CCC1'), ('Cyclopentane', 'C1CCCC1'), ('Cyclohexane', 'C1CCCCC1'),
        ('Benzene', 'c1ccccc1'), ('Naphthalene', 'c1ccc2ccccc2c1'), ('Anthracene', 'c1ccc2cc3ccccc3cc2c1'),
        ('Phenanthrene', 'c1ccc2c(c1)ccc1ccccc12'), ('Biphenyl', 'c1ccc(-c2ccccc2)cc1'),
        ('Decalin', 'C1CCC2CCCCC2C1'), ('Adamantane', 'C1C2CC3CC1CC(C2)C3'),
        ('Steroid skeleton', 'C1CCC2C(C1)CCC1C2CCC2CCCC21'),
    ],
    'Heterocycles': [
        ('Pyridine', 'c1ccncc1'), ('Pyrimidine', 'c1cncnc1'), ('Pyrazine', 'c1cnccn1'), ('Pyrrole', 'c1cc[nH]c1'),
        ('Furan', 'c1ccoc1'), ('Thiophene', 'c1ccsc1'), ('Imidazole', 'c1cnc[nH]1'), ('Pyrazole', 'c1cn[nH]c1'),
        ('Oxazole', 'c1cocn1'), ('Thiazole', 'c1cscn1'), ('Indole', 'c1ccc2[nH]ccc2c1'), ('Quinoline', 'c1ccc2ncccc2c1'),
        ('Isoquinoline', 'c1ccc2cnccc2c1'), ('Benzimidazole', 'c1ccc2[nH]cnc2c1'), ('Purine', 'c1ncc2[nH]cnc2n1'),
        ('Piperidine', 'C1CCNCC1'), ('Piperazine', 'C1CNCCN1'), ('Morpholine', 'C1COCCN1'), ('Pyrrolidine', 'C1CCNC1'),
        ('Tetrahydrofuran', 'C1CCOC1'), ('Tetrahydropyran', 'C1CCOCC1'), ('Coumarin', 'O=c1ccc2ccccc2o1'),
        ('Flavone', 'O=c1cc(-c2ccccc2)oc2ccccc12'), ('Porphine', 'c1cc2cc3ccc(cc4ccc(cc5ccc(cc1n2)[nH]5)n4)[nH]3'),
    ],
    'Amino acids (L)': [
        ('Glycine', 'NCC(=O)O'), ('Alanine', 'C[C@H](N)C(=O)O'), ('Valine', 'CC(C)[C@H](N)C(=O)O'),
        ('Leucine', 'CC(C)C[C@H](N)C(=O)O'), ('Isoleucine', 'CC[C@H](C)[C@H](N)C(=O)O'),
        ('Proline', 'O=C(O)[C@@H]1CCCN1'), ('Phenylalanine', 'N[C@@H](Cc1ccccc1)C(=O)O'),
        ('Tyrosine', 'N[C@@H](Cc1ccc(O)cc1)C(=O)O'), ('Tryptophan', 'N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O'),
        ('Serine', 'N[C@@H](CO)C(=O)O'), ('Threonine', 'C[C@@H](O)[C@H](N)C(=O)O'), ('Cysteine', 'N[C@@H](CS)C(=O)O'),
        ('Methionine', 'CSCC[C@H](N)C(=O)O'), ('Asparagine', 'NC(=O)C[C@H](N)C(=O)O'),
        ('Glutamine', 'NC(=O)CC[C@H](N)C(=O)O'), ('Aspartic acid', 'N[C@@H](CC(=O)O)C(=O)O'),
        ('Glutamic acid', 'N[C@@H](CCC(=O)O)C(=O)O'), ('Lysine', 'NCCCC[C@H](N)C(=O)O'),
        ('Arginine', 'N=C(N)NCCC[C@H](N)C(=O)O'), ('Histidine', 'N[C@@H](Cc1c[nH]cn1)C(=O)O'),
        ('L-DOPA', 'N[C@@H](Cc1ccc(O)c(O)c1)C(=O)O'),
    ],
    'Sugars': [
        ('β-D-Glucopyranose', 'OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O'),
        ('α-D-Glucopyranose', 'OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O'),
        ('β-D-Galactopyranose', 'OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O'),
        ('β-D-Fructofuranose', 'OC[C@H]1O[C@@](O)(CO)[C@@H](O)[C@@H]1O'),
        ('β-D-Ribofuranose', 'OC[C@H]1O[C@@H](O)[C@H](O)[C@@H]1O'),
        ('Sucrose', 'OC[C@H]1O[C@@](CO)(O[C@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)[C@@H](O)[C@@H]1O'),
    ],
    'Nucleobases': [
        ('Adenine', 'Nc1ncnc2[nH]cnc12'), ('Guanine', 'Nc1nc2[nH]cnc2c(=O)[nH]1'), ('Cytosine', 'Nc1cc[nH]c(=O)n1'),
        ('Thymine', 'Cc1c[nH]c(=O)[nH]c1=O'), ('Uracil', 'O=c1cc[nH]c(=O)[nH]1'),
    ],
    'Betalain chemistry': [
        ('Betalamic acid', 'OC(=O)C1CC(=CC=O)C=C(N1)C(=O)O'),
        ('cyclo-DOPA', 'O=C(O)[C@@H]1Cc2cc(O)c(O)cc2N1'),
        ('L-DOPA', 'N[C@@H](Cc1ccc(O)c(O)c1)C(=O)O'),
    ],
    'Solvents and reagents': [
        ('Water', 'O'), ('Methanol', 'CO'), ('Ethanol', 'CCO'), ('Acetonitrile', 'CC#N'), ('Acetone', 'CC(C)=O'),
        ('DMSO', 'CS(C)=O'), ('DMF', 'CN(C)C=O'), ('Dichloromethane', 'ClCCl'), ('Chloroform', 'ClC(Cl)Cl'),
        ('Ethyl acetate', 'CCOC(C)=O'), ('THF', 'C1CCOC1'), ('Toluene', 'Cc1ccccc1'), ('Formic acid', 'O=CO'),
        ('Acetic acid', 'CC(=O)O'), ('Trifluoroacetic acid', 'OC(=O)C(F)(F)F'), ('Triethylamine', 'CCN(CC)CC'),
        ('Pyridine', 'c1ccncc1'), ('Sulfuric acid', 'O=S(=O)(O)O'), ('Nitric acid', 'O=[N+]([O-])O'),
    ],
}
