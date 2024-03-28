"""
    c3sdb/build_utils/_parsing.py

    Dylan Ross (dylan.ross@pnnl.gov)

    module with various parsing utilities
"""


import re
from typing import Optional, Dict, List


# list of names to exclude that technically match the peptide regex but are not actually peptides
PEPTIDE_EXCLUSION_LIST: List[str] = [
    'ADP', 'AMP', 'ATP', 'CDP', 'CMP', 'CTP', 'GDP', 'GMP', 'GTP', 'UMP', 'UDP', 'UTP', 'IMP',
    'EDTA', 'DEET', 'AICAR',
    'ADMA', 'SDMA',
    'SAH', 'SAM', 'NADP', 'NADPH',
    'EDDP', 'LSD', 'MDEA', 'MDMA', 'MDPV',
    'DGDP', 'DGTP', 'FAD', 'FMN', 'ALLIIN', 'ICARIIN', 'TADALAFIL', 'MELPHALAN', 'THIRAM', 
    'VALSARTAN', 'ASPARTAME', 'KAWAIN', 'ASARYLALDEHYDE', 'EPICATECHIN', 'ANIRACETAM', 
    'PINACIDIL', 'CEFDINIR',
    'HMMA', 'PYRIMETHANIL', 'ACETAMIPRID', 'THC', 'FENTANYL', 'ENALAPRIL',
    'RISPERIDAL', 'PFNA', 'PFDA', 'PFNA', 'PFDA', 'PMPA', 'EPN', 'MCPA',
    'APIIN'
]


# regex patterns to exclude from peptide labeled compounds
PEPTIDE_EXCLUSION_PATS: List[str] = [
    r'.*I[DN]E$',
    r'.*ATE$',
    r'.*[CNTRGDLSPLV]IN$',
    r'^M*ETH.*',
    r'.*[AE]NE$',
    r'.*EIN$'
]


def parse_carbohydrate(name: str
                       ) -> bool :
    """
    Checks a compound name for matches with regex patterns for carbohydrates, returns
    a boolean reflecting whether the name looks like a carbohydrate.

    patterns: 

    - '[(]Hex[)][0-9]+'
    - '^([LD][-])*[A-Za-z]+ose$'
    - '^[A-Za-z]+itol$'
    - '.*Cyclodextrin'
    - 'Lacto[-]N[-][A-Za-z]+ose.*'
    
    Parameters
    ----------
    name : ``str``
        compound name
    
    Returns
    -------
    matches : ``bool``
        flag indicating if the name matches one of the carbohydrate
        regex patterns
    """
    if re.match(r'[(]Hex[)][0-9]+', name):
        return True
    if re.match(r'^([LD][-])*[A-Za-z]+ose$', name):
        return True
    if re.match(r'^[A-Za-z]+itol$', name):
        return True
    if re.match(r'.*Cyclodextrin', name):
        return True
    if re.match(r'Lacto[-]N[-][A-Za-z]+ose.*', name):
        return True
    return False


def parse_lipid(name: str
                ) -> Optional[Dict[str, str | int]] :
    """
    parses a lipid name into lipid class and fatty acid composition, returning a
    dictionary with the information. 
    
    .. note::

        Handles total fatty acid composition, as well
        as individual composition, examples:
            PC(38:3)        --> class: PC, n_carbon: 38, n_unsat: 3
            PC(18:1/20:2)   --> class: PC, n_carbon: 38, n_unsat: 3, 
                                fa_comp: ((n_carbon: 18, n_unsat: 1), (n_carbon: 20, n_unsat: 2))
        Also, handles special fatty acid notations (modifiers) used for ceramides and 
        plasmalogen lipids, examples:
            Cer(d36:2)      --> class: Cer, n_carbon: 36, n_unsat: 2, fa_mod: d
            Cer(d18:1/18:1) --> class: PC, n_carbon: 38, n_unsat: 3, fa_mod: d,
                                fa_comp: ((n_carbon: 18, n_unsat: 1), (n_carbon: 18, n_unsat: 1))
            PE(p40:4)       --> class: PE, n_carbon: 40, n_unsat: 4, fa_mod: p
            PE(p20:2/20:2)  --> class: PE, n_carbon: 40, n_unsat: 4, fa_mod: p,
                                fa_comp: ((n_carbon: 20, n_unsat: 2), (n_carbon: 20, n_unsat: 2))
    
    .. note:: 
    
        lipid name must conform to the general format:
        <lipid_class>([modifier]<n_carbon>:<n_unsat>[/<n_carbon>:<n_unsat>[/<n_carbon>:<n_unsat>]])
    
    Parameters
    ----------
    name : ``str``
        compound name
    
    Returns
    -------
    lipid_info : ``dict(...)`` or None
        parsed lipid information (always contains 'class', 'n_carbon', and 'n_unsat'
        attributes) or None if it cannot be parsed as a lipid
    """
    parsed = {}
    # compile regex pattern
    s = (r"^(?P<cls>[A-Za-z]+)\((?P<mod>[pdoe]*)(?P<fc1>[0-9]+):(?P<fu1>[0-9]+)"
         r"[/_]*((?P<fc2>[0-9]+):(?P<fu2>[0-9]+))*[/_]*((?P<fc3>[0-9]+):(?P<fu3>[0-9]+))*\)")
    l_pat = re.compile(s)
    # parse the name using regex
    l_res = l_pat.match(name)
    if l_res:
        # lipid class (required)
        if l_res.group('cls'):
            parsed["lipid_class"] = l_res.group('cls')
        else:
            #msg = "parse_lipid: failed to parse lipid class for: {}".format(name)
            #raise ValueError(msg)
            return None      
        # fc1 and fu1 are always required
        if not l_res.group('fc1') or not l_res.group('fu1'):
            #raise_fatty_acid_value_error()
            return None
        # check if a second fatty acid composition is supplied, e.g. (18:1/16:0)
        # if so, need to compute total fatty acid composition and add individual
        # fatty acids to a list
        if l_res.group('fc2'):
            if not l_res.group('fu2'):
                #raise_fatty_acid_value_error()
                return None
            # add info from the first two fatty acid compositions
            fc1, fu1 = int(l_res.group('fc1')), int(l_res.group('fu1'))
            fc2, fu2 = int(l_res.group('fc2')), int(l_res.group('fu2'))
            parsed["fa_comp"] = [
                {"n_carbon": fc1, "n_unsat": fu1}, 
                {"n_carbon": fc2, "n_unsat": fu2}
            ]
            # check for 3rd FA composition
            fc3, fu3 = 0, 0
            if l_res.group('fc3'):
                if not l_res.group('fu3'):
                    #raise_fatty_acid_value_error()
                    return None
                fc3, fu3 = int(l_res.group('fc3')), int(l_res.group('fu3'))
                parsed["fa_comp"].append({"n_carbon": fc3, "n_unsat": fu3})
            # compute total fatty acid composition
            parsed["n_carbon"] = fc1 + fc2 + fc3
            parsed["n_unsat"] = fu1 + fu2 + fc3
        else:
            # fc1 and fu1 are the total fatty acid composition
            parsed["n_carbon"] = int(l_res.group('fc1'))
            parsed["n_unsat"] = int(l_res.group('fu1'))

        # add fatty acid modifier if present
        if l_res.group('mod'):
            parsed["fa_mod"] = l_res.group('mod')
    else:
        # could not parse name as a lipid
        parsed = None
    return parsed


def parse_peptide(name: str, 
                  exclude: bool = True
                  ) -> bool :
    """
    Checks a compound name for matches with a regex for peptide sequences,
    optionally excluding names that look like peptides but are not really
    
    Parameters
    ----------
    name : ``str``
        compound name
    exclude : ``bool``, default=True
        exclude names from a statically defined peptide exclusion list 
        (`PEPTIDE_EXCLUSION_LIST`) or names that match any defined regex patterns
        to exclude (`PEPTIDE_EXCLUSION_PATS`)
    
    Returns
    -------
    matches : ``bool``
        flag indicating if the name matches the peptide regex pattern
    """
    expats = [
        re.compile(pat) for pat in PEPTIDE_EXCLUSION_PATS
    ]    
    if exclude:
        excluded = False
        for expat in expats:
            if expat.match(name):
                excluded = True
        excluded = name in PEPTIDE_EXCLUSION_LIST or excluded
    pep_pat = re.compile(r'^[ACDEFGHIKLMNPQRSTVWY]+$')
    pep_res = pep_pat.match(name)
    if pep_res and not excluded:
        return True
    return False

