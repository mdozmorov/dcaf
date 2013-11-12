"""
General statistical procedures. 

Most of these procedures accept and return
:py:class:`pandas.DataFrame` or :py:class:`pandas.Series` objects, or
their sparse equivalents. The term "matrix" in this context means
a :py:class:`pandas.DataFrame`.
"""

def standardize(X):
    """
    Standardize a matrix such that each column has a mean of 0 
    and standard deviation of 1.

    :param X: The matrix
    :type X: :py:class:`pandas.DataFrame`
    :returns: The standardized matrix
    :rtype: :py:class:`pandas.DataFrame`
    """
    return X.apply(lambda x: (x - x.mean()) / x.std(), axis=1)
