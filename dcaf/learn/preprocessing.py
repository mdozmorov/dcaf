"""
Transformers and algorithms for preprocessing data for machine learning.

Most of the classes here use the scikits-learn fit/transform API.
"""

import numpy

__all__ = ["InvalidColumnRemover", "QuantileNormalizer"]

class SimpleTransform(object):
    def fit_transform(self, X, y=None):
        self.fit(X,y)
        return self.transform(X,y)

class InvalidColumnRemover(SimpleTransform):
    """
    A transformer that removes columns with too many missing values.
    """
    def __init__(self, min_count=5):
        assert(min_count >= 0)
        self._min_count = min_count
    
    def fit(self, X, y=None):
        self._valid_columns = numpy.invert(numpy.isnan(X)).sum(axis=0) \
                              >= self._min_count
        return self

    def transform(self, X, y=None):
        return X[:,self._valid_columns]
        
class QuantileNormalizer(SimpleTransform):
    """
    A transformer that quantile normalizes samples. 
    
    It assumes there are no missing values. If this is not the case,
    an imputation or NaN removal step needs to be placed earlier in the pipeline.
    
    For more information about quantile normalization, particularly in the
    context of expression analysis, see:

    - http://bioinformatics.oxfordjournals.org/content/19/2/185
    """
    __slots__ = ["_copy", "coef_"]

    def __init__(self, copy=True):
        self._copy = copy

    def fit(self, X, y=None):
        self.coef_ = X.mean(axis=0)
        self.coef_.sort()
        return self
        
    def transform(self, X, y=None):
        if self._copy:
            X = X.copy()
        for i in range(X.shape[0]):
            ix = X[i,:].argsort()
            X[i,ix] = self.coef_
        return X
