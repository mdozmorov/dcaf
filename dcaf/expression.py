"""
Manipulate transcript expression data.
"""
import argparse
import sys
import numpy
import pandas

from sklearn.preprocessing import Imputer

import dcaf.util

from dcaf.db import DCAFConnection
from dcaf.io import read_matrix

standardize = lambda X: X.apply(lambda x: (x - x.mean()) / x.std(), axis=1)

def pearson_distance(X, Y):

    # Calculate the distance matrix
    n = len(X.columns) - 1
    D = 1 - (numpy.dot(standardize(X), standardize(Y).T) / n)
    return pandas.DataFrame(D, index=X.index, columns=Y.index)

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
        X = prepare(X.T)
    else:
        X = prepare(read_matrix(args.x))

    Y = prepare(read_matrix(args.y) if (args.x and args.y) else X)

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
        D.to_csv(sys.stdout, float_format="%0.3f", sep="\t")
    else:
        D[numpy.isnan(D)] = numpy.inf
        for i,rn in enumerate(X.index):
            print(rn,*(Y.index[j] for j in \
                    D.iloc[i,:].argsort()[:args.k]), sep="\t")

if __name__ == "__main__":
    pairwise_distance()
