# Development Notes

##### __2021-04-13__

_Updates:_
*  Added perfluoro- compounds dataset from Dodds _et al._ (src_tag: `dodd0220`), there are only about 50 compounds but they are worth including because there are not many perfluoro- compounds in the database.
* Ran database setup, there are now __14056__ compounds in the database, __13065__ (92.9 %) of which have MQNs.

---
##### __2019-09-16__

_Updates:_
*  Added some more of the edited SMILES structures, mostly some more weird bacterial lipids
* Ran database setup, there are now __6492__ compounds in the database, __6250__ (96.3 %) of which have SMILES structures.

---
##### __2019-09-16__

_Updates:_
*  Cleaned the E. Coli Lipid data from Kelly and Libin's Chemistry and Physics of Lipids Paper (src_tag: `hine0119`). 179 CCS values for bacterial lipids, all with manually generated SMILES structures. 
* tried running build script, had an issue with improperly escaped backslashes in the SMILES structures. Fixed it.
* Ran database setup, there are now __6492__ compounds in the database, __6221__ (95.8 %) of which have SMILES structures. 

---
##### __2019-08-31__

_Updates:_
*  Libin and Jang Ho did some work to fill in SMILES structures still missing from the database. They also made some edits to some bacterial lipids incorrectly labeled as LPG when they are actually LysPG (_i.e._, lysyl-PG). I am beginning to apply some of those changes to the database.
* Step 1: added SMILES structures to the mislabeled bacterial lipids (in the cleaned .json data: `hine1217.json`) and change the labels accordingly
* Step 2: Some compounds were identified that were mislabeled with chemical class. Added a build script (`fix_mislabeled_chem_class.py`) to be run directly after `label_chem_class.py` that fixes specific instances of mislabeled compounds. There are also several suggestions for names to change, that may need to be a separate script that runs earlier. Not planning on calling this script in the main build script until it is implemented fully... 

---
##### __2019-07-20__

_Updates:_
* Started analysis on C3S.db ... performed PCA, PLS-RA. Did ML trials on the data ... these are the primary results for CCSML paper

---
##### __2019-07-17__

_Updates:_
* we need a column in the master table to store mass (in addition to m/z that we already have) in order to handle compounds that have multiple charges (CCS vs. mass plots have more consistent trends than CCS vs. m/z plots when multiple charges are present). While I am at it, I will also add a charge (z) column. The charge can easily be obtained from the adduct information (_e.g._ `[M+3H]3+` obviosly has a charge of 3) and mass is just m/z * z.
* added `mass` and `z` columns to `C3SDB_schema.sql` and adjusted the `fill_db_from_src.py` build script to compute these values from the `adduct` and `mz` columns and add them in when building the database. 
* added a file `db_setup/src_pdfs/reference_info.py` containing data structures for connecting `src_tag` with actual reference information, complete with links to the publication websites.
* added a table to the database (defined in `CCS_pred_schema.sql`) for representing predicted CCS values and associated information. Upon training a predictive model, predicted CCS values can be added to this table for the reference data (those data that are contained in this database) and for future user-supplied data. Since we will be using some untargeted classification on the front end of the CCS prediction, there is a column for the class label as well (this will just be integers, since the classification task is untargeted). Adding this table will occur at the end of the build script (`build_new_db.sh`)
* built the (rough) chemical class labels into the database. Added a new built script to do this (`label_chem_class.py`) which is executed after MQNs are generated. Added a `carbohydrate_parser.py` utility script that does the same thing that `peptide_parser.py` does but for carbohydrates. Made a bunch of changes to `peptide_parser.py` so that certain compound names that look like they could be peptides will not be improperly parsed as peptides.
* re-ran the build scripts with the classification in place:

`chem_class_label` | count (%)
------- | -------
small molecule | 4705 (77.0 %)
lipid | 1329 (21.7 %)
carbohydrate | 167 (2.7 %)
peptide | 112 (1.8 %)

---
##### __2019-06-23__

_Updates:_
* finished docstrings for functions in `C3SData.db_interface` module
* turned db_setup into a module (added an `__init__.py`) to make documentation generation easier and to link all of the individual build script documentation linked together conveniently
* set up documentation build script and directory for `db_setup` module, generated documentation for everything
* made a copy of `db_setup/generate_lipid_smiles.py` (`db_setup/generate_lipid_smiles_v2.py`) with a much better implementation. By grouping similar lipid classes into a small number of bigger functions, I was able to reduce a ton of repetative code. Also implemented lysophospholipids and some glycerolipids while I was at it.
* completely replaced `db_setup/generate_lipid_smiles.py` with `db_setup/generate_lipid_smiles_v2.py` since I got it to work and it had better class coverage anyways
* re-ran DB setup, there are now __6313__ entries in the `master` table of `C3S.db`, __6042__ (95.7 %) of which have SMILES structures 


---
##### __2019-06-12__

_Updates:_
* added the carbohydrate and TAA data from the J. C. May, et al. 2014 dataset along with manually obtained SMILES structures. 
* There are now __6313__ entries in the `master` table of `C3S.db`, __5935__ of which have SMILES structures 

---
##### __2019-06-09__

_Updates:_
* added another dataset (Stow, _et al._, _Anal. Chem._ __2017__) with 80-something drift tube values of peptides and other metabolomics types of compounds in it. Source tag is `stow0817`.
* there was something wrong with the SMILES structures generated for peptides preventing production of MQNs, sequences were being terminated with 'OH' but I think it should only be O (implicit hydrogens unless considering adducts or stereochemistry)... made that change and tried to build the database again... That worked! Now only a handfull of compounds having SMILES structures did not get MQNs (SMILES with errors, all from `bijl0517`)... fixed!

_C3S.db stats:_

group | count (percent)
------- | -------
all compounds  | 6204 (100 %)
`smi` NOT NULL | 5826 (94 %)
missing `smi`  | 378 (6 %) 
have MQNs      | 5826 (94 %)
missing MQNs   | 378 (6 %)

_missing `smi` stats:_

`src_tag` | count (percent missing `smi`)
------- | -------
`hine1217` | 106 (28 %)
`zhen0917` | 80 (21 %)
`may_0114` | 58 (15 %)
`nich1118` | 54 (14 %)
`righ0218` | 28 (7 %)
`hine0217` | 28 (7 %)
`zhou1016` | 12 (3 %)
`zhou0817` | 7 (2 %)
`stow0817` | 5 (1 %)

* The main source of missing SMILES structures is from the `hine1217` dataset. These are bacterial lipids so they are not well represented in the LIPID MAPS database and I do not have any functions implemented for generating generic SMILES structures.
* Added columns to the `master` table, `ccs_type` and `ccs_method`, reflecting the nature of the measurements in the combined collection. This includes information like "DT" or "TW" for the type of mobility measurement and other method information like calibrants.

_C3S.db stats:_

group | count (percent)
------- | -------
all compounds  | 6204 (100 %)
drift tube measurements | 3594 (58 %)
traveling-wave measurements | 2610 (42 %)

--- 
##### __2019-04-23__

_Updates:_
* added a `gen_docs.sh` script to the `C3SData/doc` directory, uses `pydoc` to generate HTML documentation from docstrings then moves everything into the `C3SData/doc` directory
* finished descriptive docstrings for `C3SD.featurize(...)`, `C3SD.train_test_split(...)`, and `C3SD.center_and_scale(...)`
* implemented `C3SD.center_and_scale(...)`
* added and implemented `C3SD.get_categorical_y(...)` method, needed for `C3SD.train_test_split(...)` to stratify on rough CCS distribution. Uses the following binning scheme based on the full CCS distribution (quartiles):
```
                Q1            Q2             Q3
                |             |              |
          bin1  | bin2 | bin3 | bin4 | bin 5 | bin 6
                       |             |
              (Q2 - Q1 / 2) + Q1     |
                            (Q3 - Q2 / 2) + Q2
```
* implemented `C3SD.train_test_split(...)`, seems to be working for both stratify on CCS and dataset

--- 
##### __2019-04-09__

_Updates:_
* Created a python package, `C3SData`, for interfacing with the data in `C3S.db` to provide easy access, filtering, train/test set splitting, etc. for use in building/evaluating ML models for CCS prediction. 
* Began defining the main class (`C3SData.data.C3SD`) that will actually be used to access data. Wrote detailed docstrings defining the target behavior for the `C3SD.__init__(...)`, `C3SD.featurize(...)`, `C3SD.train_test_split(...)`, and `C3SD.center_and_scale(...)` methods.
* added a module (`C3SData/db_interface.py`) for doing all of the SQL work/actually accessing the database. Added functions `fetch_single_dset(...)` and `fetch_multi_dset(...)` and started work on their docstrings.
* added a `docs` directory, no documentation generation scripts yet...
* Began implementing `C3SD.__init__(...)`, `fetch_single_dset(...)`, and `fetch_multi_dset(...)`. So far all seem to be working as expected now
* added `generate_mqns.py` to the build script and made sure it was working. There are a number of SMILES structures that were unable to be interpreted by `rdkit` and probably just need minor fixes like removing random spaces and stuff like that. 
* Implemented `C3SD.featurize(...)`, seems to be working properly.

--- 
##### __2019-04-08__

`build_new_db.sh` currently calls the following build scripts:
* `fill_db_from_src.py` 
* `fetch_smiles.py`

_Updates:_
* Created `peptide_parser.py` script to use in the `fetch_lipid_smiles.py` workflow. Uses a simple regex (`'^[ACDEFGHIKLMNPQRSTVWY]+$'`) to match names that look like peptides. 
* Created `generate_peptide_smiles.py` script to generate SMILES structures from linear peptide sequences.
* Fixed an issue in `generate_lipid_smiles.py` where lipid classes 'SM' and 'Cer' were being ignored if they had the fatty acid modifier 'd', which is the wrong behavior because that modifier is already accounted for in the 'SM' and 'Cer' SMILES generator functions.
* Assorted fixes to chemical names in `src_zhou1016.json` to produce some more matches when searching PubChem.
* Lots of assorted fixes to chemical names in `src_zhen0917.json` to produce some more matches when searching PubChem/LipidMAPS and to make lipid SMILES generation work better

After running `fetch_smiles.py` __381__ of the 6118 total entries (6.2%) are still missing SMILES structures:

`src_tag` | entries without SMILES
------- | -------
`hine0217`|28
`hine1217`|106
`may_0114`|58
`nich1118`|58
`righ0218`|28
`zhen0917`|84
`zhou0817`|7
`zhou1016`|12

---

##### __2019-01-17__

`build_new_db.sh` currently calls the following build scripts:
* `fill_db_from_src.py` 
* `fetch_smiles.py`

Made some modifications to the `generate_lipid_smiles.py` build script:
* Added functions for generating generic SMILES structures for PCs and PEs. 
* Modified the `generate_lipid_smiles` function so that when trying to match the class, the fatty acid modifier is taken into account (`if lipid_cls in lipid_cls_smi_func:` → `if lipid_cls in lipid_cls_smi_func and not fa_mod:`). So far fatty acid modifiers are left unimplemented, but I think it would be best to add an `fa_mod` kwarg to the function calls for relevant lipid class SMILES generator functions and just pass it along rather than ignoring entries with fatty acid modifications entirely. 


After running `fetch_smiles.py` __579__ of the 6118 total entries (9.5%) are still missing SMILES structures:

`src_tag` | entries without SMILES
------- | -------
`hine0217` | 62
`hine1217` | 106
`may_0114` | 147
`nich1118` | 63
`righ0218` | 30
`zhen0917` | 136
`zhou0817` | 7
`zhou1016` | 28

The number of entries missing SMILES structures has actually _gone up_ with these fixes, which indicates that we were erroneously generating SMILES structures for lipids containing fatty acid modifiers as though the modifiers were not present. All the more reason to actually implement fatty acid modifiers in the individual generator functions.

The overall structure of `generate_lipid_smiles.py` is pretty messy and has a bunch of repeated code. I think it would be better to implement some sort of base function(s) or decorators to make them more modular. The only real differences from function to function is the `base_smi` and the manner in which the `nc` and `nu` are calculated from the fatty acid composition. The whole script could use some refactoring.  

--- 
##### __2019-01-16__

`build_new_db.sh` currently calls the following build scripts:
* `fill_db_from_src.py` 
* `fetch_smiles.py`

After running `fill_db_from_src.py` the database contains __6118__ entries:

`src_tag` | entries
------- | -------
`bijl0517` | 205
`groe0815` | 132
`hine0217` | 257
`hine0817` | 1426
`hine1217` | 163
`may_0114` | 389
`moll0218` | 357
`nich1118` | 1092
`pagl0314` | 96
`righ0218` | 106
`zhen0917` | 597
`zhou0817` | 451
`zhou1016` | 847

After running `fetch_smiles.py` __545__ of the 6118 total entries (8.9%) are still missing SMILES structures:

`src_tag` | entries without SMILES
------- | -------
`hine0217` | 14
`hine0817` | 4
`hine1217` | 106
`may_0114` | 157
`nich1118` | 63
`righ0218` | 30
`zhen0917` | 136
`zhou0817` | 7
`zhou1016` | 28

Most of these entries that are still missing SMILES structures are lipids for which a search of LIPID MAPS returned no results and a generic SMILES structure could not be generated (mostly due to unimplemented lipid classes in `generate_lipid_smiles.py`). Another significant source of compounds missing SMILES structures are the peptide data from the `may_0114` dataset. I need to figure out a good way of generating peptide SMILES structures from amino acid sequence. I fixed the 4 compounds from the `hin0817` dataset by manually adding entries for benzalkonium chlorides (C16, C14, C12) to `search_cache.json` and changing the compound name "TOREMIPHENE CITRATE" → "TOREMIFENE CITRATE" so that it should be recognized when searching against PubChem. 

---
