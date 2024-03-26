"""
    kmcm.py
    Dylan H. Ross
    2019/08/02

    description:
        TODO
"""


from sklearn.base import BaseEstimator, RegressorMixin, clone
from sklearn.cluster import KMeans
from numpy import argwhere, array
from itertools import product



class KMCMulti(BaseEstimator, RegressorMixin):
    """
KMCMulti
    description:
        TODO
"""
    
    def __init__(self, seed=69, n_clusters=3, use_estimator=None, estimator_params=None):
        """
KMCMulti.__init__
    description:
        TODO
    parameters:
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

    def fit(self, X, y):
        """
KMCMulti.fit
    description:
        TODO
    parameters:
        X (array-like) -- features 
        y (array-like) -- targets
    returns: 
        () -- 
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
        
    def predict(self, X):
        """
KMCMulti.predict
    description:
        TODO
    parameters:
        X (array-like) -- features
    returns:
        (numpy.ndarray) -- predictions 
"""
        return array([self.estimators_[self.kmeans_.predict(x.reshape(1, -1))[0]].predict(x.reshape(1, -1))[0] for x in X])

    
def kmcm_p_grid(n_clusters, est_params):
    """
kmcm_p_grid
    description:
        generates a parameter grid that can be used with GridSearchCV for hyperparameter tuning
    parameters:
        n_clusters (list(int)) -- list of values to try for n_clusters
        est_params (dict(str:list(...))) -- values to try for individual estimator parameters, in
                                            the style of the parameter grid used for GridSearchCV
    returns:
        (list(dict(...))) -- parameter grid for use with GridSearchCV
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
    

