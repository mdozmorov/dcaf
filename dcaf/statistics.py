"""
General statistical procedures. 

Most of these procedures accept and return
:py:class:`pandas.DataFrame` or :py:class:`pandas.Series` objects, or
their sparse equivalents. The term "matrix" in this context means
a :py:class:`pandas.DataFrame`.
"""

def scale(X):
    """
    Scale a matrix such that each column has a mean of 0 
    and standard deviation of 1. 
    
    See also:

    - `:py:func:sklearn.preprocessing.scale`

    :param X: The matrix
    :type X: :py:class:`pandas.DataFrame`
    :returns: The standardized matrix
    :rtype: :py:class:`pandas.DataFrame`
    """
    return (X - X.mean(axis=0, skipna=True)) / X.std(axis=0, skipna=True)
