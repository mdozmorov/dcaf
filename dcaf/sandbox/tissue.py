import os

import pandas
import numpy
import pandas
import sklearn.decomposition

from sklearn.cross_validation import KFold, cross_val_score
from sklearn.preprocessing import StandardScaler, Imputer
from sklearn.pipeline import Pipeline
from sklearn.metrics import log_loss
from sklearn.multiclass import OneVsRestClassifier
from sklearn.svm import SVC

from dcaf.db import DCAFConnection
from dcaf.expression import impute_expression
from dcaf.statistics import standardize
from dcaf.util import coo_to_df

URSA_MEAN = 7.32
URSA_STDEV = 2.63

def tissue_expression_training_set(taxon_id=9606, limit=100):
    db = DCAFConnection.from_configuration()
    c = db("""
    SELECT sample_term.sample_id, expression.data, 
        sample_term.term_id, sample_term.probability
    FROM sample_term
    INNER JOIN term
    ON term.id=sample_term.term_id
    INNER JOIN ontology
    ON ontology.id=term.ontology_id
    INNER JOIN sample
    ON sample.id=sample_term.sample_id
    INNER JOIN expression
    ON expression.sample_id=sample.id
    INNER JOIN platform
    ON sample.platform_id=platform.id
    INNER JOIN taxon
    ON platform.taxon_id=taxon.id
    WHERE ontology.namespace='BTO'
    AND sample_term.probability=1
    AND taxon.id=%s
    ORDER BY random()
    LIMIT %s""", taxon_id, limit)
    samples, data, tissues, values = zip(*c)
    T = coo_to_df(zip(samples, tissues, values))
    T.index.name = "Sample ID"
    T.columns.name = "Term ID"
    c.execute("""SELECT id FROM gene WHERE gene.taxon_id=%s ORDER BY id""", 
              (taxon_id,))
    X = pandas.DataFrame.from_records(list(data),
                                      index=samples, 
                                      columns=[e[0] for e in c])
    return X,T


def ursa_infer_tissue(X):
    """
    Infer tissue for sample(s) using URSA. The columns should be gene symbols.
    """
    if X.columns.dtype != 'O':
        raise Exception("The expression matrix columns are not object types, implying they do not contain gene symbols. Gene symbols are required for tissue imputation.")

    URSA_BIN = which("ursa")
    Xt = impute(X).T
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "24"
    with tempfile.NamedTemporaryFile() as h:
        X.T.dropna().to_csv(h.name, sep="\t")
        p = Popen([URSA_BIN, "-i", h.name, "-n", str(X.shape[0]), "-r", "1"],
                  cwd=dirname(URSA_BIN), stdout=PIPE, env=env)
        T = pandas.io.parsers.read_csv(p.stdout, skiprows=(6*X.shape[0] + 1),
                                       sep="\t", index_col=0, header=False,
                                       names=X.index)
    c = db("""
    SELECT translate(term.name, ' ', '_'), term.id
    FROM term
    INNER JOIN ontology
    ON ontology.id=term.ontology_id
    WHERE ontology.namespace='BTO'""")
    m = dict([(k.replace("/", "--"),v) for k,v in c])
    T.index = [m[k] for k in T.index]
    T.index.name = "Term ID"
    T.columns.name = "Sample ID"
    return T

# TODO: implement generic disk-based memoization
def load_training_set():
    if not os.path.exists("scratch/X.pkl"):
        X, T = tissue_expression_training_set(limit=100000)
        X.to_pickle("scratch/X.pkl")
        T.to_pickle("scratch/T.pkl")
    else:
        X = pandas.read_pickle("scratch/X.pkl")
        T = pandas.read_pickle("scratch/T.pkl")
    return X,T

def default_multilabel_pipeline():
    pca = sklearn.decomposition.PCA(5)
    inner = SVC(probability=True)
    ovr = OneVsRestClassifier(inner)
    model = Pipeline(steps=[
        ('scale', StandardScaler()),
        ('impute', Imputer(axis=0)),
        ('pca', pca), 
        ('ovr', ovr)
    ])
    return model

class MultilabelClassifier(object):
    """
    A convenience wrapper for multilabel classification allowing
    input and output of `:py:class:pandas.DataFrame` objects.
    """
    def __init__(self, pipeline=default_multilabel_pipeline()):
        self._pipeline = pipeline

    def fit(self, X, y):
        pass

def multilabel_cv(X,Y):
    Y = Y.ix[:,Y.sum() > 5]
    model = multilabel_pipeline()

    # List of score strings:
    # http://scikit-learn.org/stable/modules/model_evaluation.html#model-evaluation
    score = "log_loss"
    return cross_val_score(model, X, y=Y.to_dense().fillna(0).astype(int),
                           scoring=score, verbose=True)

def bto_term_names():
    db = DCAFConnection.from_configuration()
    return dict(db("""
    SELECT term.id, term.name 
    FROM term
    INNER JOIN ontology
    ON term.ontology_id=ontology.id
    WHERE ontology.namespace='BTO'"""))

if __name__ == "__main__":
    main()
