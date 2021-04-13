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
- [ ] add a `fa_mod` kwarg to the individual lipid class SMILES generator functions in `generate_lipid_smiles.py` to be able to handle fatty acid modifiers. The initial implementation can just ignore them and return `None`, thus achieving the same behavior as in the current implementation where any lipid containing a fatty acid modifier is just ignored in the `generate_lipid_smiles` function call.
- [ ] assess whether any of the SMILES structures fetched from PubChem (or other databases) contain any extra noncovalent molecules like salts (indicated by a `.` in the SMILES string). If so, need to write a small script that fixes those so that the SMILES only contains the compound of interest. 


__C3SData__
- [ ] implement `C3SD.__repr__(...)` and/or `C3SD.__str__(...)` methods for debugging purposes

Development notes [here](dev_notes.md).
