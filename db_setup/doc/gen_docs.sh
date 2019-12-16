#!/bin/bash

python3 -m pydoc -w db_setup
python3 -m pydoc -w db_setup.fetch_smiles
python3 -m pydoc -w db_setup.fill_db_from_src
python3 -m pydoc -w db_setup.generate_lipid_smiles
python3 -m pydoc -w db_setup.generate_mqns
python3 -m pydoc -w db_setup.generate_peptide_smiles
python3 -m pydoc -w db_setup.lipid_parser
python3 -m pydoc -w db_setup.peptide_parser
python3 -m pydoc -w db_setup.carbohydrate_parser
python3 -m pydoc -w db_setup.lmaps_scrape
python3 -m pydoc -w db_setup.pubchem_scrape
python3 -m pydoc -w db_setup.label_chem_class


mv -f *.html db_setup/doc/
