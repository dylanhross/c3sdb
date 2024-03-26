"""
    peptide_parser.py
    Dylan Ross
    2019/04/08

        Utility for parsing peptide names -- simple regex match for the pattern: 
            '^[ACDEFGHIKLMNPQRSTVWY]+$' 
"""


import re


def parse_peptide(name, exclude=True):
    """
parse_peptide
    description:
        Checks a compound name for matches with a regex for peptide sequences, returns
        a boolean reflecting whether the name looks like a peptide.
    parameters:
        name (str) -- compound name
        [exclude (bool)] -- exclude matches to names in the exclusion_list [optional, default=True] 
    returns:
        (bool) -- name matches peptide regex
"""
    # list of names to exclude that technically match the peptide regex but are not actually peptides
    exclusion_list = [
        'ADP', 'AMP', 'ATP', 'CDP', 'CMP', 'CTP', 'GDP', 'GMP', 'GTP', 'UMP', 'UDP', 'UTP', 'IMP',
        'EDTA', 'DEET', 'AICAR',
        'ADMA', 'SDMA',
        'SAH', 'SAM', 'NADP', 'NADPH',
        'EDDP', 'LSD', 'MDEA', 'MDMA', 'MDPV',
        'DGDP', 'DGTP', 'FAD', 'FMN', 'ALLIIN', 'ICARIIN', 'TADALAFIL', 'MELPHALAN', 'THIRAM', 'VALSARTAN',
        'ASPARTAME', 'KAWAIN', 'ASARYLALDEHYDE', 'EPICATECHIN', 'ANIRACETAM', 'PINACIDIL', 'CEFDINIR'
    ]

    expat1 = re.compile(r'.*I[DN]E$')
    expat2 = re.compile(r'.*ATE$')
    expat3 = re.compile(r'.*[CNTRGDLSPLV]IN$')
    expat4 = re.compile(r'^M*ETH.*')
    expat5 = re.compile(r'.*[AE]NE$')
    expat6 = re.compile(r'.*EIN$')

    if exclude:
        excluded = name in exclusion_list or expat1.match(name) or expat2.match(name) or expat3.match(name) \
                     or expat4.match(name) or expat5.match(name) or expat6.match(name)
    else:
        excluded = False

    pep_pat = re.compile(r'^[ACDEFGHIKLMNPQRSTVWY]+$')
    pep_res = pep_pat.match(name)
    if pep_res and not excluded:
        return True
    return False
