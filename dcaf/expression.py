"""
Manipulate transcript expression data.

In this module, a "matrix" or "expression matrix" is a
:py:class:`pandas.DataFrame` which contains transcripts as columns and
experiments as rows. 
""" 

import argparse 
import sys 
import numpy
import pandas

from sklearn.preprocessing import Imputer

import dcaf.util
import dcaf.statistics

from dcaf.db import DCAFConnection
from dcaf.io import read_matrix

def impute_expression(X):
    """
    Impute the value of missing genes to be the mean expression of that gene
    across all experiments.
    
    Also drops experiments and genes with no data.

    :param X: The expression matrix, containing NAs for missing values
    :type X: :py:class:`pandas.DataFrame`
    :returns: An expression matrix with NAs replaced by imputed values, 
      and with rows or columns with insufficient data for imputation removed
    :rtype: :py:class:`pandas.DataFrame`
    """
    if not pandas.isnull(X).sum().sum():
        return X
    X = X.dropna(axis=(0,1), how="all")
    model = Imputer(axis=1, strategy="mean")
    return DataFrame(model.fit_transform(X),
                     index=X.index,
                     columns=X.columns)

def pearson_distance(X, Y=None): 
    """
    Return a generator containing, for each row in X, the Pearson
    distance between that row and each row in Y, if Y is provided.
    Otherwise, the self-distance between each element in X.
    
    :param X: The first matrix
    :type X: :py:class:`pandas.DataFrame`
    :param Y: The second matrix
    :type Y: :py:class:`pandas.DataFrame`
    :returns: A generator yielding an element for each row in 
      X describing the Pearson distance between that row and each row in Y 
    :rtype: generator of :py:class:`pandas.Series`
    """
    if Y is None:
        Y = X
    assert(numpy.array(X.columns == Y.columns).all())
    n = len(X.columns) - 1
    X = standardize(X)
    Y = standardize(Y)
    for i in X.index:
        dist = 1 - numpy.dot(Y, X.ix[i,:])
        yield pandas.Series(dist, index=Y.index)

@dcaf.util.entry_point
def pairwise_distance(argv):
    """
    Given one or two matrices, calculate a distance metric
    between each pair of rows, and output either the k nearest
    neighbors for each element, or the entire distance matrix.

    The default output format is as follows:
    - one tab-delimited line per row in the first ("X") matrix
    - the first element in the line is the "X" row ID
    - the remaining N elements (modifiable by the -n switch) \
      are the N rows in "Y" which have the smallest distance \
      to the current row in "X"

    Alternatively, the -m flag will allow the complete distance matrix
    to be output.

    The default metric is Pearson distance. 
    """
    # FIXME: Implement other distance metrics. 
    #   See: http://gedas.bizhat.com/dist.htm
    # FIXME: Implement -n for files.

    parser = argparse.ArgumentParser(
            description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-x", type=argparse.FileType("r"),
            help="Path to the first matrix.")
    parser.add_argument("-y", type=argparse.FileType("r"),
            help="""Path to the second matrix (optional). 
            If this is not provided, X will be compared against itself.""")
    parser.add_argument("-t", "--taxon-id", type=int,
                        help="""Instead of using matrix files, 
                        query the DCAF database for an expression matrix
                        from this taxon ID.""")

    parser.add_argument("-k", type=int, default=50,
        help="Number of top correlated transcripts to print.")
    parser.add_argument("-n", "--n-samples", type=int, default=5000,
                        help="""The maximum number of (randomly selected) 
                        samples to use for similarity calculation. 
                        Can be useful when RAM is limited.""")

    parser.add_argument("--output_matrix", '-m', action="store_true",
            default=False,
            help="Output the entire distance matrix.")

    args = parser.parse_args(argv)


    def prepare(M):
        M = M.dropna(how='all', axis=(0,1))
        return pandas.DataFrame(Imputer().fit_transform(M.as_matrix()),
                                index=M.index, columns=M.columns)

    if args.taxon_id:
        db = DCAFConnection.from_configuration()
        X,_ = db.expression(taxon_id=args.taxon_id, limit=args.n_samples, shuffle=True)
        X = impute_expression(X).T
    else:
        X = impute_expression(read_matrix(args.x)).T

    Y = impute_expression(read_matrix(args.y) if (args.x and args.y) else X)

    if not numpy.array(X.columns == Y.columns).all():
        raise Exception("All column IDs must match.")
    else:
        shared_cols = list(set(X.columns) & set(Y.columns))
        if not shared_cols:
            raise Exception("No columns are shared between the matrices!")
        X = X.ix[:,shared_cols]
        Y = Y.ix[:,shared_cols]

    # Calculate the distance
    D = pearson_distance(X, Y)

    # Output the results
    if args.output_matrix:
        raise NotImplementedError
        #D.to_csv(sys.stdout, float_format="%0.3f", sep="\t")
    else:
        #D[numpy.isnan(D)] = numpy.inf
        for i,d in enumerate(D):
            rn = X.index[i]
            d = d.fillna(numpy.inf)
            print(rn,*(Y.index[j] for j in d.argsort()[:args.k]) , sep="\t")

if __name__ == "__main__":
    pairwise_distance()
