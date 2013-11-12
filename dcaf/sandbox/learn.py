"""
Machine learning algorithms and helpers.
"""

class SparseMultilabelLR(object):
    def __init__(self):
        pass

class SimpleRegressor(object):
    def __init__(self, method="mean"):
        assert(method in ("mean", "mode"))
        self._method = method

    def fit(self, X, y):
        if self._method == "mean":
            self._y = y.mean()
        elif self._method == "mode":
            self._y = np.argmax(np.bincount(np.round(y).astype(int)))
    
    def predict(self, X):
        return array([self._y for _ in range(X.shape[0])])

def train(X, Y):
    """
    Train a multilabel classifier.
    - X is a DataFrame of instances (rows) versus features (columns).
    - Y is a SparseDataFrame of instances (rows) versus term IDs (columns).
    """
    #X = impute(X)
    #common = list(set(X.index) & set(Y.index))
    #X = X.ix[common,:]
    #Y = Y.ix[common,:]
    Yt = Y.T
    Yl = []
    for sample in Yt.columns:
        Yl.append(tuple(Yt.ix[:,sample].to_dense().notnull().nonzero()[0]))
    #inner = SVC(kernel="linear", probability=True)
    inner = LogisticRegression()
    #inner = SimplePrior()
    clf = OneVsRestClassifier(inner)#, n_jobs=-1)
    clf.fit(X, Yl)
    return clf

def evaluate(X, model, annotations):
    loss = []
    for gene_id in X.columns:
        M = annotations.get(int(gene_id), set())
        y = numpy.array([1 if label in M else 0 for label in model.classes_])
        y_hat = [[p, 1-p] for p in model.predict_proba(X[gene_id])[0]]
        loss.append(log_loss(y, y_hat))
    return numpy.array(loss)

from sklearn.base import BaseEstimator, ClassifierMixin
class SimplePrior(BaseEstimator, ClassifierMixin):
    "A very simple classifier that always predicts the prior probability."
    def __init__(self):
        pass
    
    def fit(self, X, Y):
        p = Y.sum() / Y.size
        self.p = numpy.array([p, 1-p])
        
    def predict_proba(self, X):
        return numpy.array([self.p,] * X.shape[0])
 
# Multilabel learning

class BRRDT(object):
    """
    Implements 
    https://www.siam.org/proceedings/datamining/2010/dm10_068_zhangx.pdf
    """
    def __init__(self, n=10):
        self._trees = []

    def fit(self, X, Y):
        pass

import tempfile
import os
import itertools
from collections import defaultdict
from subprocess import Popen, PIPE, check_output
from os.path import realpath, dirname

import numpy
from pandas import DataFrame, Series, SparseDataFrame
import pandas
from sklearn.preprocessing import Imputer
from sklearn.multiclass import OneVsRestClassifier
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import log_loss
from sklearn.linear_model import LinearRegression, Ridge, Lasso
from sklearn.svm import SVR
from numpy import array, fabs
import numpy as np
from sklearn.datasets import make_multilabel_classification
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import LabelBinarizer, Imputer
from sklearn.multiclass import OneVsRestClassifier
from sklearn.svm import LinearSVC, SVC
from sklearn.pipeline import Pipeline
from sklearn.feature_selection import SelectKBest
from pandas import DataFrame, SparseDataFrame

DataFrame.show = lambda self: self.head(n=3).T.head(n=3).T
SparseDataFrame.show = DataFrame.show

def insert_sample_terms(T, min_p=0.5):
    """
    Insert sample-term associations.
    """
    assert(T.index.name=="Term ID")
    assert(T.columns.name=="Sample ID")
    #for sample_id in T.index:
    #    c.execute("""
    #    DELETE FROM sample_term 
    #    INNER JOIN sample 
    #    ON sample.id=sample_term.sample_id
    #    WHERE sample.id=%s""", (sample_id,))
    c.executemany("""
    INSERT INTO sample_term (sample_id, term_id, probability)
    VALUES (%s,%s,%s)""", 
                  ((int(r[1]), int(r[2]), float(r[3])) for r in
                   T[T>min_p].unstack().reset_index().dropna().to_records()))
    db.commit()

def tissue_expression_training_set(taxon_id=9606, limit=200):
    c.execute("""
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
    LIMIT %s""", (taxon_id, limit))
    samples, data, tissues, values = zip(*c)
    T = coo_to_df(zip(samples, tissues, values))
    T.index.name = "Sample ID"
    T.columns.name = "Term ID"
    c.execute("""SELECT id FROM gene WHERE gene.taxon_id=%s ORDER BY id""", 
              (taxon_id,))
    X = DataFrame.from_records(list(data),
                               index=samples, columns=[e[0] for e in c])
    return X,T

       
def impute_age():
    X, P = gfa.platform_expression("GPL96")
    model = impute.KNNImputer()
    Xi = model.fit_transform(X, axis=1)

    age = array(P["age"].tolist())
    Xm = Xi.as_matrix()
    ix = array((age >= 10) & (age <= 120)).nonzero()[0]
    np.random.shuffle(ix)
    Xm = Xm[ix,:]
    age = age[ix]

    n_train = 2000
    n_test = 500
    #clf = SVR(C=1e-5, epsilon=1)
    #clf = LinearRegression()
    clf = Ridge()
    #clf = SimpleRegressor()
    #clf = Lasso()
    clf.fit(Xm[:n_train,:], age[:n_train])
    y = age[n_train:(n_train+n_test)]
    y_hat = clf.predict(Xm[n_train:(n_train+n_test)])
    dy = y - y_hat

    bias_tr = y_hat.mean() - age.mean()
    print("\nBias (vs train):\t\t", bias_tr)
    print("Bias (vs test):\t\t\t", dy.mean())
    print("Mean error:\t\t\t", fabs(dy).mean())
    print("Mean error (bias corrected):\t", fabs(dy - bias_tr).mean())
    print("MSE:\t\t\t\t", np.power(dy,2).mean())
    #clf.fit(Xi.as_matrix(), age)
 
"""
g = gfa.ontology_graph()
X,T = gfa.tissue_expression_training_set(limit=10000)
X = X.dropna(axis=1, how="all")
#from grtk.impute import KNNImputer
#Xi = KNNImputer().fit_transform(X)
Xi = DataFrame(Imputer().fit_transform(X.as_matrix()),
               index=X.index, columns=X.columns)
Xs = Xi.as_matrix()

Tt = T.T
data = {}
for i in Tt.columns:
    data[i] = {}
    for j in Tt.index[Tt.ix[:,i].to_dense().notnull()]:
        data[i][j] = 1
        if j in g:
            for ancestor in g.successors(j):
                data[i][ancestor] = 1

Ts = SparseDataFrame(data).T.to_sparse(fill_value=0)
Ts = Ts.fillna(-1).as_matrix()

Xs = Xi.as_matrix()

base_clf = Pipeline([('feature_selection', SelectKBest(k=100)),
                     ('svm', SVC(probability=True))])
clf = OneVsRestClassifier(base_clf, n_jobs=-1)
N_train = 5000
ix = np.invert((Ts[:N_train,:] == -1).all(axis=0))
clf.fit(Xs[:N_train,:], Ts[:N_train,ix])

N_test = 1000
print(roc_auc_score(Ts[N_train:(N_train+N_test),ix].flat, 
                    clf.predict_proba(Xs[N_train:(N_train+N_test),:]).flat))
"""

"""
X, P = taxon_expression(9606)
X = X.dropna(axis=1, how="all")
Xi = DataFrame(Imputer().fit_transform(X.as_matrix()),
               index=X.index, columns=X.columns)
"""

