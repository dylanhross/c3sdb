#!/usr/bin/python3


from sqlite3 import connect

from reference_info import src_tags_ordered, reference_links_ordered


con = connect('C3S.db')
cur = con.cursor()


qry = """ SELECT g_id, name, adduct, mz, mass, z, ccs, ccs_type, ccs_method, smi, src_tag FROM master """


# map src_tag to a corresponding reference link
refs = {}
for src_tag, reference_link in zip(src_tags_ordered, reference_links_ordered):
    refs[src_tag] = reference_link


data = []
for g_id, name, adduct, mz, mass, z, ccs, ccs_type, ccs_method, smi, src_tag in cur.execute(qry).fetchall():
    mz, mass, ccs = float(mz), float(mass), float(ccs)
    z = int(z)
    ccs_method = ccs_method.replace(',', '')
    smi = '' if smi is None else smi
    ref = refs[src_tag]
    data.append([g_id,name, adduct, mz, mass, z, ccs, ccs_type, ccs_method, smi, ref])


fstr = '"{}","{}","{}",{:.4f},{:.4f},{},{:.2f},"{}","{}","{}","{}"\n'
with open('ccsbase_pubchem.csv', 'w') as out:
    out.write('g_id,name,adduct,mz,mass,z,ccs,ccs_type,ccs_method,smi,ref\n')
    for cmpd in data:
        out.write(fstr.format(*cmpd))


fstr = """"{}","{}","{}","CCS_Type: '{}' | Reference: '{}'"\n"""
with open('ccsbase_pubchem_substance.csv', 'w') as out:
    out.write('PUBCHEM_EXT_DATASOURCE_REGID,PUBCHEM_SUBSTANCE_SYNONYM,PUBCHEM_EXT_DATASOURCE_SMILES,PUBCHEM_SUBSTANCE_COMMENT\n')
    for cmpd in data:
        # only need a few of the fields
        g_id,name, adduct, mz, mass, z, ccs, ccs_type, ccs_method, smi, ref = cmpd
        out.write(fstr.format(g_id, name, smi, ccs_type, ref))


