"""
    generate_lipid_smiles_v2.py
    Dylan Ross
    2019/06/23

        Utility to systematically generate SMILES structures for lipids as defined by lipid
        class and sum fatty acid composition

        * This second version is intended to improve upon my original implementation, principally by making the 
        individual lipid class SMILES generator functions more modular *
"""


def carbon_chain(c, u):
    """
carbon_chain
    description:
        Returns a SMILES representation of a carbon chain with c carbons and 
        u unsaturations
    parameters:
        c (int) -- number of carbons
        u (int) -- number of unsaturations
    returns:
        (str) -- SMILES representation of the carbon chain
"""
    smi = "C"
    c -= 1
    while u > 0:
        smi += "C=CC"
        c -= 3
        u -= 1
    while c > 0:
        smi += "C"
        c -= 1
    return smi


def sphingo_smiles(c, n_carbon, n_unsat, fa_mod):
    """
sphingo_smiles
    description:
        Returns a SMILES structure (str) for a sphingolipid using a base SMILES structure (includes head group)
    parameters:
        c (str) -- lipid class
        n_carbon (int) -- fatty acid carbons
        n_unsat (int) -- fatty acid unsaturation
        fa_mod (str) -- fatty acid modifier 
    returns:
        (str) -- SMILES structure, None if anything goes wrong
"""
    # check fatty acid modifier (should only be 'd' for this class)
    if fa_mod and fa_mod != 'd':
        # return None for unimplemented fa_mod
        return None 
    # base SMILES structure (includes head group and positions for fatty acids)
    base_smi_c = {
        'SM': '[C@](COP(=O)([O-])OCC[N+](C)(C)C)([H])(NC({})=O)[C@]([H])(O){}',
        'Cer': '[C@](CO)([H])(NC({})=O)[C@]([H])(O){}',
        'GlcCer': '[C@](CO[C@@H]1O[C@H](CO)[C@@H](O)C(O)C1O)([H])(NC({})=O)[C@]([H])(O){}',
        'HexCer': 'C(COC1OC(CO)C(O)C(O)C1O)([H])(NC({})=O)[C@]([H])(O){}'
    }
    if c not in base_smi_c:
        # the specifiec class is improper for this sphingolipid smiles generator
        e = 'sphingo_smiles: the specified class ({}) is improper for this SMILES generator'
        raise ValueError(e.format(c))
    base_smi = base_smi_c[c]
    if n_unsat == 0:
        return base_smi.format(carbon_chain(n_carbon - 19, 0), carbon_chain(15, 0))
    elif n_unsat == 1:
        return base_smi.format(carbon_chain(n_carbon - 19, 0), carbon_chain(15, 1))
    else:
        return base_smi.format(carbon_chain(n_carbon - 19, n_unsat - 1), carbon_chain(15, 1))


def phospho_smiles(c, n_carbon, n_unsat, fa_mod):
    """
phospho_smiles
    description:
        Returns a SMILES structure (str) for a phospholipid using a base SMILES structure (includes head group)
    parameters:
        c (str) -- lipid class
        n_carbon (int) -- fatty acid carbons
        n_unsat (int) -- fatty acid unsaturation
        fa_mod (str) -- fatty acid modifier 
    returns:
        (str) -- SMILES structure, None if anything goes wrong
"""
    # check fatty acid modifier 
    if fa_mod:
        # return None for unimplemented fa_mod
        return None 
    # base SMILES structure (includes head group and positions for fatty acids)
    base_smi_c = {
        'PC': 'C[N+](C)(C)CCOP(OC[C@]([H])(OC({})=O)COC({})=O)(=O)O',
        'PE': 'NCCOP(OC[C@]([H])(OC({})=O)COC({})=O)(=O)O',
        'PS': 'C(O)(=O)[C@@]([H])(N)COP(OC[C@]([H])(OC({})=O)COC({})=O)(=O)O',
        'PA': 'OP(OC[C@]([H])(OC({})=O)COC({})=O)(=O)O',
        'PG': 'OCC(O)COP(OC[C@]([H])(OC({})=O)COC({})=O)(=O)O'
    }
    if c not in base_smi_c:
        # the specifiec class is improper for this smiles generator
        e = 'phospho_smiles: the specified class ({}) is improper for this SMILES generator'
        raise ValueError(e.format(c))
    base_smi = base_smi_c[c]
    # split the fatty acid carbons and unsaturations evenly between the two fatty acids
    nc, nu = n_carbon - 2, n_unsat
    nc_a, nc_b = nc / 2, nc / 2 + nc % 2
    nu_a, nu_b = nu / 2, nu / 2 + nu % 2 
    return base_smi.format(carbon_chain(nc_a, nu_a), carbon_chain(nc_b, nu_b))


def lysophospho_smiles(c, n_carbon, n_unsat, fa_mod):
    """
lysophospho_smiles
    description:
        Returns a SMILES structure (str) for a lysophospholipid using a base SMILES structure (includes head group)
        * the generated SMILES structure corresponds to the 2-lysophospholipid *
    parameters:
        c (str) -- lipid class
        n_carbon (int) -- fatty acid carbons
        n_unsat (int) -- fatty acid unsaturation
        fa_mod (str) -- fatty acid modifier 
    returns:
        (str) -- SMILES structure, None if anything goes wrong
"""
    # check fatty acid modifier 
    if fa_mod:
        # return None for unimplemented fa_mod
        return None 
    # base SMILES structure (includes head group and positions for fatty acids)
    base_smi_c = {
        'LPC': 'C[N+](C)(C)CCOP(OC[C@]([H])(O)COC({})=O)(=O)O',
        'LPE': 'NCCOP(OC[C@]([H])(O)COC({})=O)(=O)O',
        'LPS': 'C(O)(=O)[C@@]([H])(N)COP(OC[C@]([H])(O)COC({})=O)(=O)O',
        'LPA': 'OP(OC[C@]([H])(O)COC({})=O)(=O)O',
        'LPG': 'OCC(O)COP(OC[C@]([H])(O)COC({})=O)(=O)O'
    }
    if c not in base_smi_c:
        # the specifiec class is improper for this smiles generator
        e = 'lysophospho_smiles: the specified class ({}) is improper for this SMILES generator'
        raise ValueError(e.format(c))
    base_smi = base_smi_c[c]
    # only one fatty acid position
    nc, nu = n_carbon - 1, n_unsat
    return base_smi.format(carbon_chain(nc, nu))


def glycero_smiles(c, n_carbon, n_unsat, fa_mod):
    """
glycero_smiles
    description:
        Returns a SMILES structure (str) for a glycerolipid using a base SMILES structure (includes head group)
    parameters:
        c (str) -- lipid class
        n_carbon (int) -- fatty acid carbons
        n_unsat (int) -- fatty acid unsaturation
        fa_mod (str) -- fatty acid modifier 
    returns:
        (str) -- SMILES structure, None if anything goes wrong
"""
    # check fatty acid modifier 
    if fa_mod:
        # return None for unimplemented fa_mod
        return None 
    # base SMILES structure (includes head group and positions for fatty acids)
    base_smi_c = {
        'DG': 'C(O)[C@]([H])(OC({})=O)COC({})=O',
        'DGDG': 'C(O[C@@H]1O[C@H](CO[C@H]2O[C@H](CO)[C@H](O)C(O)C2O)[C@H](O)C(O)C1O)[C@]([H])(OC({})=O)COC({})=O'
    }
    if c not in base_smi_c:
        # the specifiec class is improper for this smiles generator
        e = 'glycero_smiles: the specified class ({}) is improper for this SMILES generator'
        raise ValueError(e.format(c))
    base_smi = base_smi_c[c]
    # split the fatty acid carbons and unsaturations evenly between the two fatty acids
    nc, nu = n_carbon - 2, n_unsat
    nc_a, nc_b = nc / 2, nc / 2 + nc % 2
    nu_a, nu_b = nu / 2, nu / 2 + nu % 2 
    return base_smi.format(carbon_chain(nc_a, nu_a), carbon_chain(nc_b, nu_b))


def tg_smiles(n_carbon, n_unsat, fa_mod):
    """
tg_smiles
    description:
        Returns a SMILES structure (str) for a TG using a base SMILES structure (includes head group)
    parameters:
        c (str) -- lipid class
        n_carbon (int) -- fatty acid carbons
        n_unsat (int) -- fatty acid unsaturation
        fa_mod (str) -- fatty acid modifier
    returns:
        (str) -- SMILES structure, None if anything goes wrong
"""
    # check fatty acid modifier
    if fa_mod:
        # return None for unimplemented fa_mod
        return None
    # base SMILES structure (includes head group and positions for fatty acids)
    base_smi = '{}C(=O)OCC(COC(=O){})OC(=O){}'
    # split the fatty acid carbons and unsaturations evenly between the three fatty acids
    nc, nu = n_carbon - 3, n_unsat
    # NC needs to be at least 3 (so at least 6 in input)
    if nc < 3:
        return None
    nc_a = nc // 3 + nc % 3
    nc_b = nc // 3
    nc_c = nc // 3
    nu_a = nu // 3 + nu % 3
    nu_b = nu // 3
    nu_c = nu // 3
    return base_smi.format(carbon_chain(nc_a, nu_a), carbon_chain(nc_b, nu_b), carbon_chain(nc_c, nu_c))


def generate_lipid_smiles(lipid_cls, n_carbon, n_unsat, fa_mod=None):
    """
generate_lipid_smiles
    description:
        Generates a SMILES structure given a lipid class, sum fatty acid composition, and
        optional fatty acid modifier. Returns None if there are any errors.
    parameters:
        lipid_cls (str) -- lipid class
        n_carbon (int) -- fatty acid carbons
        n_unsat (int) -- fatty acid unsaturation
        [fa_mod (str)] -- fatty acid modifier [optional, default=None]
    returns:
        (str or None) -- SMILES structures (or None if any errors)
"""
    sphingo_cls = ['SM', 'Cer', 'GlcCer', 'HexCer']
    phospho_cls = ['PC', 'PE', 'PS', 'PA', 'PG']
    lysophospho_cls = ['LPC', 'LPE', 'LPS', 'LPA', 'LPG']
    glycero_cls = ['DG', 'DGDG']

    smi = None

    # sphingolipids
    if lipid_cls in sphingo_cls:
        smi = sphingo_smiles(lipid_cls, n_carbon, n_unsat, fa_mod)

    # phospholipids
    elif lipid_cls in phospho_cls:
        smi = phospho_smiles(lipid_cls, n_carbon, n_unsat, fa_mod)

    # lysophospholipids
    elif lipid_cls in lysophospho_cls:
        smi = lysophospho_smiles(lipid_cls, n_carbon, n_unsat, fa_mod)

    # glycerolipids
    elif lipid_cls in glycero_cls:
        smi = glycero_smiles(lipid_cls, n_carbon, n_unsat, fa_mod)

    # TG
    elif lipid_cls == 'TG':
        smi = tg_smiles(n_carbon, n_unsat, fa_mod)

    return smi


