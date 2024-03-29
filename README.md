# C3SDB (Combined Collision Cross Section DataBase)
__Dylan H. Ross__
 
- If used in published work, please cite: https://doi.org/10.1021/acs.analchem.9b05772
- CCSbase website: https://www.ccsbase.net
- Contact: Libin Xu (libinxu@uw.edu)


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

> A set of pre-trained instances of encoder, scaler and KMCM-SVR CCS prediction model
> as pickle files, and the summary figure for the prediction metrics are available in 
> the `pretrained/` directory of this repository

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
The following examples demonstrate how to train a model using the K-Means clustering with SVM
approach that was used in the original paper.

> This example assumes the database (`C3S.db`) has already been built according to directions above.

#### Example with Hyperparameter Optimization 
```python
import pickle

from sklearn.model_selection import GridSearchCV
from sklearn.svm import SVR

from c3sdb.ml.data import C3SD
from c3sdb.ml.kmcm import kmcm_p_grid, KMCMulti
from c3sdb.ml.metrics import compute_metrics_train_test, train_test_summary_figure

# intialize the dataset
data = C3SD("C3S.db", seed=2345)
data.assemble_features()
data.train_test_split("ccs")
data.center_and_scale()
# generates c3sdb_OHEncoder.pkl and c3sdb_SScaler.pkl
data.save_encoder_and_scaler()  

# this generates a ton of parameter combinations, especially with higher numbers of clusters
kmcm_svr_p_grid = kmcm_p_grid([4, 5], {"C": [1000, 10000], "gamma": [0.001, 0.1]}) 
kmcm_svr_gs = GridSearchCV(KMCMulti(seed=2345, use_estimator=SVR(cache_size=1024, tol=1e-3)),
                           param_grid=kmcm_svr_p_grid, n_jobs=-1, cv=3, scoring="neg_mean_squared_error",
                           verbose=3)
kmcm_svr_gs.fit(data.X_train_ss_, data.y_train_)
kmcm_svr_best = kmcm_svr_gs.best_estimator_

# save the trained model for use later
with open("c3sdb_kmcm_svr.pkl", "wb") as pf:
    pickle.dump(kmcm_svr_best, pf)

# compute metrics for the trained model
y_pred_train = kmcm_svr.predict(data.X_train_ss_)
y_pred_test = kmcm_svr.predict(data.X_test_ss_)
summary = compute_metrics_train_test(data.y_train_, data.y_test_, y_pred_train, y_pred_test)
train_test_summary_figure(summary, "metrics.png")
```

#### Example without Hyperparameter Optimization
```python
# ... snip ...

# this generates a ton of parameter combinations, especially with higher numbers of clusters
kmcm_svr = KMCMulti(n_clusters=5,
                    seed=2345, 
                    use_estimator=SVR(cache_size=1024, tol=1e-3), 
                    estimator_params=[
                        {"C": 10000, "gamma": 0.001} for _ in range(5)
                    ]
                    ).fit(data.X_train_ss_, data.y_train_)

# save the trained model for use later
with open("c3sdb_kmcm_svr.pkl", "wb") as pf:
    pickle.dump(kmcm_svr, pf)

# ... snip ...
```

### Inference with Trained Model
> This example assumes that the model has already been trained as described in the examples above 
> and the files `c3sdb_OHEncoder.pkl`, `c3sdb_SScaler.pkl`, and `c3sdb_kmcm_svr.pkl` are all available.
```python
import pickle

from c3sdb.ml.data import data_for_inference

# generate input data (m/z + encoded adduct + MQNs, centered and scaled) for inference
# mzs, adducts, and smis are all numpy arrays with same length
# included is a boolean mask with same shape as input arrays, indicating which rows
# are included in the generated dataset (some might fail to compute MQNs from SMILES)
X, included = data_for_inference(mzs, adducts, smis,
                                 "c3sdb_OHEncoder.pkl", "c3sdb_SScaler.pkl")

# load the trained model
with open("c3sdb_kmcm_svr.pkl", "rb") as pf:
    kmcm_svr = pickle.load(pf)

# do inference
y_pred = kmcm_svr.predict(X)
```

