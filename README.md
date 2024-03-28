# C3SDB (Combined Collision Cross Section DataBase)
__Dylan H. Ross__
 

## Overview
__C3SDB__ is an aggregation of large collision cross section (CCS) collections 
covering many diverse chemical classes and measured on a number of different 
instruments. The basic criteria for inclusion in the database are:
* sufficiently large collection of values (~100 minimum)
* publicly available
* performed in Nitrogen drift gas
* reliable measurements; general agreement between datasets on overlapping values


The purpose of this database is to be able to use large-scale CCS collections, 
either individually or in user-defined combinations, for building predictive CCS 
models using machine learning or other techniques. 


## `c3sdb` Python package
The `c3sdb` Python package has utilities for building the __C3SDB__ database from 
source datasets and training the CCS prediction model

### Setup
- This code is tested to work with `Python3.12`, I cannot guarantee compatibility with
    other Python versions
- I suggest setting up a virtual environment for installing dependencies
- Download the source code: `git clone https://github.com/dylanhross/c3sdb.git`
- `cd` into repository directory
- Install requirements: `pip3 install -r requirements.txt`
- Optional: add the repository directory to `PYTHON_PATH` environment variable so 
    that the code can be imported from the `c3sdb` package from anywhere in the filesystem 

### Building Database
- `cd` into the repository directory or make sure the repository directory is included in
    the `PYTHON_PATH` evironment variable as mentioned above
- Invoke the built-in standard database build script: `python3 -m c3sdb.build_utils.standard_build`
- This will build the database file (`C3S.db`) in the current working directory
- For more control over the build process, make a copy of the standard build script 
    (`c3sdb/build_util/standard_build.py`) in the current working directory (_e.g._ `./custom_build.py`) 
    then modify and invoke (`python3 ./custom_build.py`) that to customize the build.

### Training Prediction Model
__#TODO__
