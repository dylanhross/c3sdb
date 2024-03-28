"""
    C3SData/data.py
    Dylan H. Ross
    2019/04/09
    
    description:
        TODO
"""


from typing import List, Any, Optional
from sqlite3 import connect

from sklearn.preprocessing import LabelEncoder, OneHotEncoder, StandardScaler
from sklearn.model_selection import StratifiedShuffleSplit
from numpy import concatenate, percentile, digitize, array


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
        self.LEncoder_ = None
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

    def featurize(self, 
                  encoded_adduct: bool = True, 
                  mqn_indices: Optional[List[int]] = None
                  ) -> None :
        """
        Converts SMILES structures (self.smi_) into features for ML using a combination of m/z, encoded MS adduct 
        (using 1-hot encoding), and MQNs. The full concatenated feature set is stored in the self.X_ instance variable
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
        mqn_indices : ``List[int]`` or ``None``, default=None
            individually specify indices of MQNs to include in the feature set,
            or None to include all
        """
        ohe_adducts = None
        if encoded_adduct:
            # first encode the adducts as integers, then convert those to OneHot vectors
            # To reduce the number of adducts that have to get OneHot encoded, filter through the adducts list and
            # convert any adduct that is not among the top common adducts to a single label: 'other' 
            self.LEncoder_ = LabelEncoder()
            le_adducts = self.LEncoder_.fit_transform(self.filter_common_adducts(self.adduct_))
            self.OHEncoder_ = OneHotEncoder(sparse=False, categories='auto')
            ohe_adducts = self.OHEncoder_.fit_transform(le_adducts.reshape(-1, 1)).T
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
                dset.featurize(encoded_adduct=encoded_adduct, mqn_indices=mqn_indices)

    def filter_common_adducts(self, adducts):
        """
C3SD.filter_common_adducts
    description:
        To reduce the number of adducts that have to get OneHot encoded, filter through the adducts list and convert 
        any adduct that is not among the top common adducts to a single label: 'other' 
    parameters:
        adducts (np.ndarray(str)) -- numpy array containing adducts for the dataset
    returns:
        (np.ndarrray(str)) -- numpy array containing adducts with uncommon adducts replaced with 'other'
"""
        # all of these adducts had >100 examples in the combined CCS database
        ref = ['[M+NH4]+', '[M+Na-2H]-', '[M+K]+', '[M-H]-', '[M+Na]+', '[M+H]+']
        common = adducts.copy()
        for i in range(common.shape[0]):
            if common[i] not in ref:
                common[i] = 'other'
        return common

    def train_test_split(self, stratify, test_frac=0.2):
        """
C3SD.train_test_split
    description:
        Shuffles the data then splits it into a training set and a test set, storing each in self.X_train_,
        self.y_train_, self.X_test_, self.y_test_ instance variables. The splitting is done in a stratified manner
        based on either CCS or dataset source. In the former case, the CCS distribution in the complete dataset is
        binned into a rough histogram (8 bins) and the train/test sets are split such that they each contain similar 
        proportions of this roughly binned CCS distribution. In the latter case, the train/test sets are split such
        that they each preserve the rough proportions of all dataset sources present in the complete dataset.
        This method DOES NOT get called on C3SD objects in the self.datasets_ instance variable.

        Sets the following instance variables:
            self.X_train_       (training set split of features)
            self.y_train_       (training set split of labels)
            self.N_train_       (training set size) 
            self.X_test_        (test set split of features)
            self.y_test_        (test set split of labels)
            self.N_test_        (test set size)
            self.SSSplit_       (StratifiedShuffleSplit instance)

        ! self.featurize(...) must be called first to generate the features and labels (self.X_, self.y_) !
    parameters:
        stratify (str) -- specifies the method of stratification to use when splitting the train/test sets: 'source' 
                          for stratification on data source, or 'ccs' for stratification on CCS 
        [test_frac (float)] -- fraction of the complete dataset to reserve as a test set, defaults to an 80 % / 20 %
                               split for the train / test sets, respectively [optional, default=0.2]
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
            y_cat = self.get_categorical_y()
        # initialize StratifiedShuffleSplit
        self.SSSplit_ = StratifiedShuffleSplit(n_splits=1, test_size=test_frac, random_state=self.seed_)
        # split and store the X and y train/test sets as instance variables
        for train_index, test_index in self.SSSplit_.split(self.X_, y_cat):
            self.X_train_, self.X_test_ = self.X_[train_index], self.X_[test_index]
            self.y_train_, self.y_test_ = self.y_[train_index], self.y_[test_index]
        # store the size of the train/test sets in instance variables
        self.N_train_ = self.X_train_.shape[0] 
        self.N_test_ = self.X_test_.shape[0] 

    def get_categorical_y(self):
        """
C3SD.get_categorical_y
    description:
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
    returns:
        (np.ndarray(int)) -- categorical (binned) label data
"""
        # get the quartiles from the label distribution
        q1, q2, q3 = percentile(self.y_, [25, 50, 75])
        # get midpoints
        mp12 = q1 + (q2 - q1) / 2.
        mp23 = q2 + (q3 - q2) / 2.
        # bin boundaries
        bounds = [q1, mp12, q2, mp23, q3]
        return digitize(self.y_, bounds)

    def center_and_scale(self):
        """
C3SD.center_and_scale
    description:
        Centers and scales the training set features such that each has an average of 0 and variance of 1. Applies
        this transformation to the training and testing features, storing the results in the self.X_train_ss_ and 
        self.X_test_ss_ instance variables, respectively. Also stores a reference to the fitted StandardScaler
        instance for use with future data.
        This method DOES NOT get called on C3SD objects in the self.datasets_ instance variable.

        sets the following instance variables:
            self.SScaler_       (StandardScaler instance)
            self.X_train_ss_    (centered/scaled training set features)
            self.X_test_ss_     (centered/scaled test set features)

        ! self.train_test_split(...) must be called first to generate the training features and labels (self.X_train_, 
        self.y_train_) that are used to initialize the StandardScaler !
"""
        if self.X_train_ is None:
            msg = 'C3SD: center_and_scale: self.X_train_ is not initialized, self.train_test_split(...) must be ' + \
                  'called before calling self.center_and_scale(...)'
            raise RuntimeError(msg)

        # perform the scaling
        self.SScaler_ = StandardScaler()
        self.X_train_ss_ = self.SScaler_.fit_transform(self.X_train_)
        self.X_test_ss_ = self.SScaler_.transform(self.X_test_)







