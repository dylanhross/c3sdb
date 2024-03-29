"""
    c3sdb/ml/kmcm.py

    Dylan Ross (dylan.ross@pnnl.gov)

    Module housing the K-means clustering with multiple estimators (KMCMulti)
    ensemble estimator
"""


from itertools import product

from sklearn.base import BaseEstimator, RegressorMixin, clone
from sklearn.cluster import KMeans
from numpy import argwhere, array


class KMCMulti(BaseEstimator, RegressorMixin):
    """
    TODO: description
    """
    
    # TODO: type annotations
    def __init__(self, n_clusters, seed=69, use_estimator=None, estimator_params=None):
        """
        TODO: description

        Parameters
        n_clusters : ``int``
            number of clusters to fit
        
        TODO: finish reformatting these parameter descriptions
        [seed (int)] -- pRNG seed [optional, default=69]
        [n_clusters (int)] -- the number of clusters to fit [optional, default=3]
        [use_estimator (sklearn Regressor)] -- instance of individual estimator to use on each cluster 
                                               [optional, default=None]
        [estimator_params (list(dict(...)))] -- parameters to initialize each estimator with 
                                                [optional, default=None]
        """
        self.seed = seed
        self.n_clusters = n_clusters
        self.use_estimator = use_estimator
        self.estimator_params = estimator_params

    # TODO: type annotations
    def fit(self, X, y):
        """
        TODO: description

        TODO: reformat parameters and returns descriptions
        parameters:
        X (array-like) -- features 
        y (array-like) -- targets
        returns: 

        """
        # first fit the KMeans clustering model
        self.kmeans_ = KMeans(n_clusters=self.n_clusters, random_state=self.seed)
        self.kmeans_.fit(X)
        
        # split up the datasets by cluster
        cluster_X = [X[argwhere(self.kmeans_.labels_ == i).ravel()] for i in range(self.n_clusters)]
        cluster_y = [y[argwhere(self.kmeans_.labels_ == i).ravel()] for i in range(self.n_clusters)]
        
        # store the number of samples in each cluster
        self.cluster_sizes_ = [_.shape[0] for _ in cluster_X]
        
        # initialize individual estimators with their associated parameters
        self.estimators_ = [clone(self.use_estimator) for _ in range(self.n_clusters)]
        for est, p in zip(self.estimators_, self.estimator_params):
            est.set_params(**p)
        #self.estimators_ = [self.use_estimator(p) for _, p in zip(range(self.n_clusters), self.estimator_params)]
        
        # fit individual estimators with cluster data
        for est, cx, cy in zip(self.estimators_, cluster_X, cluster_y):
            est.fit(cx, cy)
        
        # return the fitted regressor
        return self

    # TODO: type annotations        
    def predict(self, X):
        """
        TODO: description

        TODO: finish reformatting these parameter descriptions
        Parameters
        ----------
        X (array-like) -- features
        returns:
        (numpy.ndarray) -- predictions 

        Returns
        -------
        TODO: ?
        """
        return array([self.estimators_[self.kmeans_.predict(x.reshape(1, -1))[0]].predict(x.reshape(1, -1))[0] for x in X])


# TODO: type annotations    
def kmcm_p_grid(n_clusters, est_params):
    """
    generates a parameter grid that can be used with GridSearchCV for hyperparameter tuning
    
    Parameters
    ----------
    n_clusters : ``list(int)`` -- list of values to try for n_clusters
    est_params : ``dict(str:list(...))``
        values to try for individual estimator parameters, in the style of the 
        parameter grid used for GridSearchCV
    Returns
    -------
    p_grid : ``list(dict(...))``
        parameter grid for use with GridSearchCV
"""
    # all permutations of the estimator parameters
    perms = []
    keys, values = zip(*est_params.items())
    for v in product(*values):
        perms.append(dict(zip(keys, v)))
    # parameter grid
    pg = []
    for nc in n_clusters: 
        n_perms = [perms for _ in range(nc)]
        pg.append({'n_clusters': [nc], 'estimator_params': [list(_) for _ in product(*n_perms)]})             
    return pg
    