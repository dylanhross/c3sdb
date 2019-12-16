"""
    generate_peptide_smiles.py
    Dylan Ross
    2018/10/19

        Utility to systematically generate SMILES structures for peptides using simple linear sequences
"""


# SMILES strings for each amino acid (excludes C-terminal OH)
aa_base = 'N[C@@H]({})C(=O)'
aa_to_smi = {
    'A': aa_base.format('C'),
    'C': aa_base.format('CS'),
    'D': aa_base.format('CC(=O)O'),
    'E': aa_base.format('CCC(=O)O'),
    'F': aa_base.format('CC{0}=CC=CC=C{0}'),
    'G': 'NCC(=O)',
    'H': aa_base.format('CC{0}=CNC=N{0}'),
    'I': aa_base.format('C(C)CC'),
    'K': aa_base.format('CCCCN'),
    'L': aa_base.format('CC(C)C'),
    'M': aa_base.format('CCSC'),
    'N': aa_base.format('CC(=O)N'),
    'P': 'N{0}C[C@@H](CCC{0})C(=O)',
    'Q': aa_base.format('CCC(=O)N'),
    'R': aa_base.format('CCCNC(=N)N'),
    'S': aa_base.format('CO'),
    'T': aa_base.format('C[C@@H](C)O'),
    'V': aa_base.format('C(C)C'),
    'W': aa_base.format('CC{0}=CNC{1}=C(C=CC=C{1}){0}'),
    'Y': aa_base.format('CC{0}=CC=C(O)C=C{0}')
}


def seq_to_smiles(seq):
    """
seq_to_smiles
    description:
        generates a SMILES structure for a linear peptide sequence
    parameters:
        seq (str) -- peptide sequence
    returns:
        (str) -- peptide SMILES
""" 
    smi = ''
    rings = 0
    for aa in seq:
        aa_smi = aa_to_smi[aa]
        if '{1}' in aa_smi:
            rings += 1
            r0 = rings
            rings += 1
            r1 = rings
            smi += aa_smi.format(r0, r1)
        elif '{0}' in aa_smi:
            rings += 1
            r0 = rings
            smi += aa_smi.format(r0)
        else:
            smi += aa_smi
    smi += 'O'
    return smi


