"""
    C3SData/data.py
    Dylan H. Ross
    2019/04/09
    
    description:
        TODO
"""


import os
from typing import List, Any, Optional, Tuple
from sqlite3 import connect
import pickle

from sklearn.preprocessing import OneHotEncoder, StandardScaler
from sklearn.model_selection import StratifiedShuffleSplit
# TODO: replace all of the calls to these functions with np.function(...)
#       and get rid of this itemized import in favor of the more typical
#       import numpy as np
from numpy import concatenate, percentile, digitize, array
import numpy as np
from numpy import typing as npt

from c3sdb.build_utils.mqns import compute_mqns


# define adducts with sufficient representation in the database
# to justify explicitly encoding them as features
_EXPLICIT_ADDUCTS = [
    "[M+H]+", "[M+Na]+", "[M-H]-", "[M+NH4]+", "[M+K]+",
    "[M+H-H2O]+", "[M+HCOO]-", "[M+CH3COO]-", "[M+Na-2H]-"
]


def _all_dset(db_path: str
              ) -> List[str] :
    """
    Returns a list of all src_tags in the C3S.db
    
    Parameters
    ----------
    db_path : ``str``
        path to C3S.db database file
    
    Returns
    -------
    src_tags : (list(str)) -- list of src_tags in C3S.db
    """
    con = connect(db_path)
    cur = con.cursor()
    qry = 'SELECT DISTINCT src_tag FROM master'
    src_tags = [_[0] for _ in cur.execute(qry).fetchall()]
    con.close()
    return src_tags


def _fetch_single_dset(db_path: str, 
                       src_tag: str
                       ) -> Any :
    """
    fetch name, mz, adduct, ccs, src_tag, smi, and chem_class_label for a dataset using
    the provided src_tag

    Parameters
    ----------
    db_path : ``str`` 
        path to C3S.db database file
    src_tag : ``str``
        specify the source dataset
    
    Returns
    -------
    names, mzs, adducts, ccss, srcs, smis, mqns, cls_labs : ``numpy.ndarray(...)``
    """
    con = connect(db_path)
    cur = con.cursor()
    qryA = 'SELECT g_id, name, mz, adduct, ccs, src_tag, smi, chem_class_label FROM master WHERE src_tag="{}"'.format(src_tag) + \
          ' AND smi IS NOT NULL AND g_id IN (SELECT g_id FROM mqns)'
    qryB = 'SELECT mqns.* FROM mqns INNER JOIN master ON mqns.g_id=master.g_id WHERE master.src_tag="{}"'.format(src_tag)
    names, mzs, adducts, ccss, srcs, smis, cls_labs = [], [], [], [], [], [], []
    for _, name, mz, adduct, ccs, src, smi, cls_lab in cur.execute(qryA).fetchall():
        names.append(name)
        mzs.append(mz)
        adducts.append(adduct)
        ccss.append(ccs)
        srcs.append(src)
        smis.append(smi)
        cls_labs.append(cls_lab)
    # fetch mqns separately from the rest of the data
    mqns = []
    for _, *mqn in cur.execute(qryB).fetchall():
        mqns.append(mqn)
    con.close()
    # self.cmpd_, self.mz_, self.adduct_, self.ccs_, self.src_, self.smi_
    return array(names), array(mzs), array(adducts), array(ccss), array(srcs), array(smis), array(mqns), array(cls_labs)


def _fetch_multi_dset(db_path: str, 
                      src_tags: List[str]
                      ) -> Any : 
    """
    fetch name, mz, adduct, ccs, src_tag, smi, and chem_class_label for a set of datasets using
    the provided src_tag

    Parameters
    ----------
    db_path : ``str`` 
        path to C3S.db database file
    src_tags : ``list(str)``
        specify the source dataset
    
    Returns
    -------
    names, mzs, adducts, ccss, srcs, smis, mqns, cls_labs : ``numpy.ndarray(...)``
    """
    con = connect(db_path)
    cur = con.cursor()
    qryA = 'SELECT g_id, name, mz, adduct, ccs, src_tag, smi, chem_class_label FROM master WHERE src_tag IN ('
    qryB = 'SELECT mqns.* FROM mqns INNER JOIN master ON mqns.g_id=master.g_id WHERE master.src_tag IN ('
    for tag in src_tags:
        qryA += '"{}",'.format(tag)
        qryB += '"{}",'.format(tag)
    qryA = qryA.rstrip(',') + ') AND smi IS NOT NULL AND g_id IN (SELECT g_id FROM mqns)'
    qryB = qryB.rstrip(',') + ')'
    names, mzs, adducts, ccss, srcs, smis, cls_labs = [], [], [], [], [], [], []
    for _, name, mz, adduct, ccs, src, smi, cls_lab in cur.execute(qryA).fetchall():
        names.append(name)
        mzs.append(mz)
        adducts.append(adduct)
        ccss.append(ccs)
        srcs.append(src)
        smis.append(smi)
        cls_labs.append(cls_lab)
    # fetch mqns separately from the rest of the data
    mqns = []
    for _, *mqn in cur.execute(qryB).fetchall():
        mqns.append(mqn)
    con.close()
    # self.cmpd_, self.mz_, self.adduct_, self.ccs_, self.src_, self.smi_, self.mqn_
    return array(names), array(mzs), array(adducts), array(ccss), array(srcs), array(smis), array(mqns), array(cls_labs)


def _filter_common_adducts(adducts: npt.ArrayLike
                           ) -> npt.ArrayLike :
        """
        To reduce the number of adducts that have to get OneHot encoded, filter through an array of
        adducts and change any adduct that is not among common adducts (defined in _EXPLICIT_ADDUCTS
        constant) to a single label: 'other' 
        
        Parameters
        ----------
        adducts : ``np.ndarray(str)``
            numpy array containing adducts for the dataset
        
        Returns
        -------
        common_adducts : ``np.ndarrray(str)``
            numpy array containing adducts with uncommon adducts replaced with 'other'
        """
        common = adducts.copy()
        for i in range(common.shape[0]):
            if common[i] not in _EXPLICIT_ADDUCTS:
                common[i] = 'other'
        return common


class C3SD:
    """
    Object for interfacing with the C3S.db database and retrieving data from it. Responsible for filtering data
    on various criteria (notably, data source), producing features for ML, handling test/train set splitting, 
    and data transfomations (normalization/centering/scaling)
    """

    def __init__(self, 
                 db_path: str, 
                 datasets: str | List[str] = [], 
                 seed: int = 69
                 ) -> None :
        """
        Initializes a new C3SD object using the path to the C3S.db database file. Uses the datasets specified in the
        datasets parameter to filter the database and build a combined dataset. If no datasets parameter is defined, 
        the default behavior is to use all of the datasets in the database. 
        The database path and pRNG seed are stored in instance variables, respectively:
        - self.db_path_
        - self.seed_
        
        The compound names, m/z, MS adduct, CCS, dataset source, SMILES structures, MQNs, and (rough) class labels are 
        all fetched and stored as numpy.ndarray in instance variables, respectively: 
        - self.cmpd_
        - self.mz_ 
        - self.adduct_
        - self.ccs_ 
        - self.src_
        - self.smi_
        - self.mqn_
        - self.cls_lab_
        
        An instance variable is created to hold individual datasets that make up the combined dataset. In the case of
        datasets=None or datasets=[...], this will be a list of C3SD objects, each containing individual datasets. Some
        methods will propagate transformations to the objects in this list, while others will not. In the case of an 
        initialization with datasets='...', this will simply be the str provided with the datasets parameter and 
        certain methods will become unavailable.
        - self.datasets_
        
        The total number of compounds in the dataset is stored in an instance variable:
        - self.N_
        
        The following instance variables are initialized as None (must be set by calls to other methods):
        - self.X_             (full array of features -> set by self.featurize(...))
        - self.y_             (full array of labels -> set by self.featurize(...))
        - self.n_features_    (number of features -> set by self.featurize(...))
        - self.LEncoder_      (LabelEncoder instance -> set by self.featurize(...))
        - self.OHEncoder_     (OneHotEncoder instance -> set by self.featurize(...))
        - self.X_train_       (training set split of features -> set by self.train_test_split(...))
        - self.y_train_       (training set split of labels -> set by self.train_test_split(...))
        - self.N_train_       (training set size -> set by self.train_test_split(...)) 
        - self.X_test_        (test set split of features -> set by self.train_test_split(...))
        - self.y_test_        (test set split of labels -> set by self.train_test_split(...))
        - self.N_test_        (test set size -> set by self.train_test_split(...))
        - self.SSSplit_       (StratifiedShuffleSplit instance -> set by self.train_test_split(...))
        - self.SScaler_       (StandardScaler instance -> set by self.center_and_scale(...))
        - self.X_train_ss_    (centered/scaled training set features -> set by self.center_and_scale(...))
        - self.X_test_ss_     (centered/scaled test set features -> set by self.center_and_scale(...))

        Parameters
        ----------
        db_path : ``str``
            path to C3S.db database file  
        datasets ``str`` or ``list(str)``, default=[]
            a list of datasets to include in the combined dataset to be fetched 
            from the database, if empty list include all datasets, if a str then 
            fetch only a single dataset 
        seed : ``int``, default=69
            pRNG seed to use for any data preparation steps with a stochastic component, stored in the
            self.seed_ instance variable 
        """
        # store database file path and pRNG seed
        self.db_path_, self.seed_ = db_path, seed
        # fetch data from the database
        if type(datasets) == str:
            # fetch a single dataset
            self.cmpd_, self.mz_, self.adduct_, self.ccs_, self.src_, self.smi_, self.mqn_, self.cls_lab_ = _fetch_single_dset(self.db_path_, datasets) 
            self.datasets_ = datasets
        elif type(datasets) == list:
            if not datasets:
                datasets = _all_dset(self.db_path_)
            # fetch either a subset of datasets or all datasets
            self.cmpd_, self.mz_, self.adduct_, self.ccs_, self.src_, self.smi_, self.mqn_, self.cls_lab_ = _fetch_multi_dset(self.db_path_, datasets)
            self.datasets_ = [C3SD(self.db_path_, datasets=dset, seed=self.seed_) for dset in datasets]
        # total number of compounds
        self.N_ = self.cmpd_.shape[0]
        # declare instance variables to use later
        self.X_ = None
        self.y_ = None
        self.n_features_ = None
        self.OHEncoder_ = None
        self.X_train_ = None
        self.y_train_ = None
        self.N_train_ = None
        self.X_test_ = None
        self.y_test_ = None
        self.N_test_ = None
        self.SSSplit_ = None
        self.SScaler_ = None
        self.X_train_ss_ = None
        self.X_test_ss_ = None

    def assemble_features(self, 
                          encoded_adduct: bool = True, 
                          mqn_indices: Optional[str | List[int]] = "all"
                          ) -> Any :
        """
        Assembles features for ML using a combination of m/z, encoded MS adduct (using 1-hot encoding), and MQNs. 

        The full concatenated feature set is stored in the self.X_ instance variable
        and self.ccs_ is simply copied into self.y_ instance variable. The self.n_features_ instance variable is also
        set accordingly. 
        The 1-hot encoder for MS adducts is stored in the self.OHEncoder_ instance variable (or just None if 1-hot 
        encoding is being used)
        Optionally, encoded_adduct may be set to False to omit the encoded adduct from the feature set.
        As another option, an array containing the indices of specific MQNs to include may be provided, if an
        empty array is provided the MQNs are omitted entirely. 
        If encoded_adduct=False and mqn_indices=[] only the m/z is used in the feature set.
        If the self.datasets_ instance variable is a list (of C3SD objects), then featurize is also called for each of
        these objects using the same parameters. 

        Sets the following instance variables:
        - self.X_             (full array of features)
        - self.y_             (full array of labels)
        - self.n_features_    (number of features)
        - self.LEncoder_      (LabelEncoder instance)
        - self.OHEncoder_     (OneHotEncoder instance)
        
        Parameters
        ----------
        encoded_adduct : ``bool``, default=True
            include encoded MS adduct in the feature set
        mqn_indices : ``str`` or ``List[int]`` or ``None``, default="all"
            individually specify indices of MQNs to include in the feature set,
            or "all" to include all or None to exclude
        """
        ohe_adducts = None
        if encoded_adduct:
            # convert adducts to OneHot vectors
            # To reduce the number of adducts that have to get OneHot encoded, filter through the adducts list and
            # convert any adduct that is not among the top common adducts to a single label: 'other' 
            self.OHEncoder_ = OneHotEncoder(sparse_output=False, categories='auto')
            common_adducts = _filter_common_adducts(self.adduct_).reshape(-1, 1)
            ohe_adducts = self.OHEncoder_.fit_transform(common_adducts).T
        use_mqns = None
        if mqn_indices:
            if mqn_indices == 'all':
                mqn_indices = [_ for _ in range(42)]
            # index self.mqn_ to get the specified indices
            use_mqns = self.mqn_.T[mqn_indices]
        # start with m/z, then concatenate ohe_adducts (if set) and use_mqns (if set) 
        x = self.mz_.reshape(1, -1)
        if ohe_adducts is not None:
            x = concatenate([x, ohe_adducts])
        if use_mqns is not None:
            x = concatenate([x, use_mqns])
        # set self.X_ and self.y_, and n_features
        self.X_ = x.T
        self.n_features_ = self.X_.shape[1]
        self.y_ = self.ccs_
        # call featurize on all of the C3SD objects in self.datasets_ (if there are any)
        if type(self.datasets_) == list:
            for dset in self.datasets_:
                dset.assemble_features(encoded_adduct=encoded_adduct, mqn_indices=mqn_indices)

    def train_test_split(self, 
                         stratify: str, 
                         test_frac: float = 0.2
                         ) -> None :
        """
        Shuffles the data then splits it into a training set and a test set, storing each in self.X_train_,
        self.y_train_, self.X_test_, self.y_test_ instance variables. The splitting is done in a stratified manner
        based on either CCS or dataset source. In the former case, the CCS distribution in the complete dataset is
        binned into a rough histogram (8 bins) and the train/test sets are split such that they each contain similar 
        proportions of this roughly binned CCS distribution. In the latter case, the train/test sets are split such
        that they each preserve the rough proportions of all dataset sources present in the complete dataset.
        This method DOES NOT get called on C3SD objects in the self.datasets_ instance variable.

        Sets the following instance variables:
        - self.X_train_       (training set split of features)
        - self.y_train_       (training set split of labels)
        - self.N_train_       (training set size) 
        - self.X_test_        (test set split of features)
        - self.y_test_        (test set split of labels)
        - self.N_test_        (test set size)
        - self.SSSplit_       (StratifiedShuffleSplit instance)

        .. note:: 
            
            self.featurize(...) must be called first to generate the features and labels (self.X_, self.y_)

        Parameters
        ----------
        stratify : ``str``
            specifies the method of stratification to use when splitting the train/test sets: 'source' 
            for stratification on data source, or 'ccs' for stratification on CCS 
        test_frac : ``float``, default=0.2
            fraction of the complete dataset to reserve as a test set, defaults to an 80 % / 20 %
            split for the train / test sets, respectively
        """
        # make sure self.featurize(...) has been called
        if self.X_ is None:
            msg = 'C3SD: train_test_split: self.X_ is not initialized, self.featurize(...) must be called before ' + \
                    'calling self.train_test_split(...)'
            raise RuntimeError(msg)
        # make sure stratify is a valid option (if provided)
        if stratify not in ['source', 'ccs']:
            msg = 'C3SD: train_test_split: stratify="{}" invalid, must be "source" or "ccs"'.format(stratify)
            raise RuntimeError(msg)
        if stratify == 'source':
            # stratify on dataset source
            y_cat = self.src_
        else:
            y_cat = self._get_categorical_y()
        # initialize StratifiedShuffleSplit
        self.SSSplit_ = StratifiedShuffleSplit(n_splits=1, test_size=test_frac, random_state=self.seed_)
        # split and store the X and y train/test sets as instance variables
        for train_index, test_index in self.SSSplit_.split(self.X_, y_cat):
            self.X_train_, self.X_test_ = self.X_[train_index], self.X_[test_index]
            self.y_train_, self.y_test_ = self.y_[train_index], self.y_[test_index]
        # store the size of the train/test sets in instance variables
        self.N_train_ = self.X_train_.shape[0] 
        self.N_test_ = self.X_test_.shape[0] 

    def _get_categorical_y(self
                           ) -> Any :
        """
        transforms the labels into 'categorical' data (required by StratifiedShuffleSplit) by performing a binning
        operation on the continuous label data. The binning is performed using the rough distribution of label values
        with the following bounds (based on quartiles):
                Q1            Q2             Q3
                |             |              |
          bin1  | bin2 | bin3 | bin4 | bin 5 | bin 6
                       |             |
              (Q2 - Q1 / 2) + Q1     |
                            (Q3 - Q2 / 2) + Q2
        Uses labels stored in the self.y_ instance variable
        
        Returns
        -------
        cat_y : ``np.ndarray(int)``
            categorical (binned) label data
        """
        # get the quartiles from the label distribution
        q1, q2, q3 = percentile(self.y_, [25, 50, 75])
        # get midpoints
        mp12 = q1 + (q2 - q1) / 2.
        mp23 = q2 + (q3 - q2) / 2.
        # bin boundaries
        bounds = [q1, mp12, q2, mp23, q3]
        return digitize(self.y_, bounds)

    def center_and_scale(self
                         ) -> None :
        """
        Centers and scales the training set features such that each has an average of 0 and variance of 1. Applies
        this transformation to the training and testing features, storing the results in the self.X_train_ss_ and 
        self.X_test_ss_ instance variables, respectively. Also stores a reference to the fitted StandardScaler
        instance for use with future data.
        This method DOES NOT get called on C3SD objects in the self.datasets_ instance variable.

        sets the following instance variables:
            self.SScaler_       (StandardScaler instance)
            self.X_train_ss_    (centered/scaled training set features)
            self.X_test_ss_     (centered/scaled test set features)

        .. note:: 
            
            self.train_test_split(...) must be called first to generate the training features and 
            labels (self.X_train_, self.y_train_) that are used to initialize the StandardScaler
        """
        if self.X_train_ is None:
            msg = 'C3SD: center_and_scale: self.X_train_ is not initialized, self.train_test_split(...) must be ' + \
                  'called before calling self.center_and_scale(...)'
            raise RuntimeError(msg)
        # perform the scaling
        self.SScaler_ = StandardScaler()
        self.X_train_ss_ = self.SScaler_.fit_transform(self.X_train_)
        self.X_test_ss_ = self.SScaler_.transform(self.X_test_)

    def save_encoder_and_scaler(self,
                                encoder_f: str = "c3sdb_OHEncoder.pkl",
                                scaler_f: str = "c3sdb_SScaler.pkl"
                                ) -> None :
        """
        save the fitted instances of the encoder and scaler objects in pickle format
        for loading later

        Parameters
        ----------
        encoder_f : ``str``, default="c3sdb_OHEncoder.pkl"
        scaler_f : ``str``, default="c3sdb_SScaler.pkl"
            file names to save the encoder and scaler instances to, respectively
        """
        # ensure both the encoder and scaler instances are present
        # TODO: check for these and raise RuntimeErrors with descriptive error messages
        assert self.OHEncoder_ is not None
        assert self.SScaler_ is not None
        with open(encoder_f, "wb") as pf:
            pickle.dump(self.OHEncoder_, pf)
        with open(scaler_f, "wb") as pf:
            pickle.dump(self.SScaler_, pf)


def data_for_inference(mzs: npt.ArrayLike, 
                       adducts: npt.ArrayLike, 
                       smis: npt.ArrayLike, 
                       encoder_f: str, 
                       scaler_f: str
                       ) -> Tuple[npt.NDArray[np.float64], npt.NDArray[np.bool_]] :
    """
    generate data for inference using lists of m/zs, adducts, smiles structures
    and fitted instances of OneHotEncoder and StandardScaler

    .. note::
        
        This function assumes you are just using the full feature set of m/z + encoded 
        adduct + all standard 42 MQNs. This function would need to be modified or other
        functions would need to be implemented in order to support other feature sets.

    Parameters
    ----------
    mzs : ``arraylike(float)``
        input m/zs
    adducts : ``arraylike(str)``
        input adducts (will be run through `_filter_common_adducts` first to ensure 
        they are all encodable)
    smis : ``arraylike(str)``
        input SMILES structures
    encoder_f : ``str``
    scaler_f : ``str``
        paths to pickle files with
        fitted instances of OneHotEncoder and StandardScaler for encoding adducts and 
        scaling features, respectively
    
    Returns
    -------
    data : ``numpy.ndarray(float)``
        scaled data ready for inference
    included : ``numpy.ndarray(bool)``
        array of booleans (same shape as input arrays) indicating which rows from the 
        original input data were able to have MQNs computed for them and thus were 
        included in the output data array
    """
    # ensure the encoder and scaler pickle files exist
    # TODO: run these checks and raise appropriate errors with descriptive error messages
    assert os.path.isfile(encoder_f)
    assert os.path.isfile(scaler_f)
    # load the encoder and the scaler instances from pickle files
    with open(encoder_f, "rb") as pf:
        encoder = pickle.load(pf)
    with open(scaler_f, "rb") as pf:
        scaler = pickle.load(pf)
    # ensure mzs, adducts, and smis have the same shape
    # TODO: run these checks and raise ValueErrors with descriptive error messages
    assert len(mzs) == len(adducts)
    assert len(adducts) == len(smis)
    # filter and encode the adducts
    enc_adducts = encoder.transform(_filter_common_adducts(np.array(adducts)).reshape(-1, 1))
    # add features row-by-row, skip any for which generating MQNs fails
    features = []
    included = []
    for mz, enc_adduct, smi in zip(mzs, enc_adducts, smis):
        if (mqns := compute_mqns(smi)) is not None:
            features.append([mz] + enc_adduct.tolist() + mqns)
            included.append(True)
        else:
            included.append(False)
    return scaler.transform(np.array(features)), included


