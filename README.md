# C3SDB (Combined Collision Cross Section DataBase)
__Dylan H. Ross__
 

## Overview
__C3SDB__ is an aggregation of large collision cross section (CCS) collections covering many diverse chemical classes and measured on a number of different instruments. The basic criteria for inclusion in the database are:
* sufficiently large collection of values (~100 minimum)
* publicly available
* performed in Nitrogen drift gas
* reliable measurements; general agreement between datasets on overlapping values


The purpose of this database is to be able to use large-scale CCS collections, either individually or in user-defined combinations, for building predictive CCS models using machine learning or other techniques. 


## Interface and Usage
_TODO_

## TODOs

__db_setup__
- [x] expand the `generate_lipid_smiles.py` script to cover more classes (_e.g._ PC, PE, LPC)
- [ ] add a `fa_mod` kwarg to the individual lipid class SMILES generator functions in `generate_lipid_smiles.py` to be able to handle fatty acid modifiers. The initial implementation can just ignore them and return `None`, thus achieving the same behavior as in the current implementation where any lipid containing a fatty acid modifier is just ignored in the `generate_lipid_smiles` function call.
- [x] refactor individual lipid class SMILES generator functions to be cleaner and less repetitive ~~(maybe use function decoration?)~~
- [x] ~~build case insensitivity into matching compound names from the search cache?~~ it works fine as is 
- [x] find a good way of generating peptide SMILES structures from amino acid sequence ~~(preferably using some type of web API, but it sounds like `openbabel` can do this so maybe better than scraping)~~ it was actually easy enough to just implement a simple function for systematic SMILES generation
- [ ] assess whether any of the SMILES structures fetched from PubChem (or other databases) contain any extra noncovalent molecules like salts (indicated by a `.` in the SMILES string). If so, need to write a small script that fixes those so that the SMILES only contains the compound of interest. 
- [x] ~~add a table to the database for ionized SMILES structures?~~ no, just keep the neutral SMILES and adduct information together
- [x] add MQNs table to database. Computing the MQNs during database build is necessary because it moves the rdkit dependency to the build scripts and out of the `C3SData` package. Get the `generate_mqns.py` build script working.
- [x] generate documentation for the build scripts, add a `doc` directory to put them into within the `db_setup` directory
- [x] add columns to the database schema for CCS instrument/technique type (call it `ccs_type`): 'DT' for drift tube or 'TW' for traveling-wave (down the road, maybe 'TJ' for modeled trajoectory method or 'ML' for predicted CCS from machine learning). Also include a `ccs_method` column for describing how the measurement was taken (_e.g._ 'stepped-field' for a DT CCS value or 'calibrated with polyalanine' for a TW CCS value). These attributes can easily be added into the build scripts, each source will have a specific set of these descriptors that apply to their values.
- [x] The J. C. May, et al. 2014 dataset seems to be missing all of its carbohydrate (and TAA) CCS values? there only seems to be the lipid and peptide values in the cleaned data... need to add those in
- [x] need to add in some rough class labels in the `master` table of `C3S.db`. This is probably most easily done during the step that SMILES structures are assigned, since there is already logic in place to distinguish lipids or peptides from the general small molecules by name. That, or copy that set of logic from the `fetch_smiles.py` build script and make a new build script dedicated to labeling the data with classes.
- [ ] ~~Add in a secondary step/script for assigning the rough chemical classes that can go in an fix specific examples of misclassified compounds. Probably easiest to do with a simple dictionary of specific compound names and their correct class label. This might be best as a standalone script (`fix_rough_class_label.py` or something like that)~~ added in the script, the idea sounds fine but just need to actually implement it


__C3SData__
- [x] Create a module for interfacing with the C3S.db and handling data. This module will be responsible for fetching data from the database and filtering it on various criteria (importantly, by source). It will also handle feature generation and train/test set splitting.
- [x] create a module in `C3SData` that handles all of the `SQLite3` querying
- [x] write descriptive docstrings for `C3SD.featurize(...)`, `C3SD.train_test_split(...)`, and `C3SD.center_and_scale(...)` (will inform implementation)
- [x] flesh out `C3SD.featurize(...)`, `C3SD.train_test_split(...)`, and `C3SD.center_and_scale(...)` implementations as described in their docstrings
- [x] make a `docs` directory in `C3SData` and put together a script for automated documentation generation/updating using `pydoc`
- [x] docstrings for functions in `C3SData.db_interface` module
- [ ] implement `C3SD.__repr__(...)` and/or `C3SD.__str__(...)` methods for debugging purposes

__analysis__
- [x] Need to generate a model using the typical _Scikit-Learn_ API that performs K-means clustering on the raw inputs, then uses the results to partition the dataset into a small number of subsets, each of which having its own estimator used to actually get the CCS prediction. In a sense, the K-means clustering would decide first if a compound 'looks like' a lipid/protein/etc. then shuttle it along to the appropriate estimator to get a CCS value. The _Scikit-Learn_ API can be achieved by inheriting from the `BaseEstimator` and `RegressorMixin` classes and implementing the appropriate methods. Using such an API allows for easy integration into the _Scikit-Learn_ ecosystem (hyperparameter tuning via `GridSearchCV`, computing perfomance metrics, _etc._). More info [here](https://scikit-learn.org/stable/developers/contributing.html#rolling-your-own-estimator).
- [x] Analyze the composition and characteristics of the full database. Compute a PCA and see how the rough chemical class labels map onto the PCA projections (_i.e._ are chemical classes a major source of variance in the overall dataset?). Determine which MQNs are invariant in the dataset (and therefore, can be omitted). Determine the most influential MQNs with respect to the variance in the dataset as a whole, as well as with respect to variance in CCS (_i.e._ loadings in PCA, PLS-DA, respectively. Also for CCS, correlation between individual MQNs and CCS). Compare DTIM and TWIM CCS values as well, see if they show any differences.

## Development Notes

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
