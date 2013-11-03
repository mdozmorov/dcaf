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

import grtk.db

db = grtk.db.connect()
c = db.cursor()
DataFrame.show = lambda self: self.head(n=3).T.head(n=3).T
SparseDataFrame.show = DataFrame.show

def which(exe):
    return realpath(check_output(["which", exe]).strip())

def multimap(pairs):
    m = {}
    for k,v in pairs:
        if not k in m:
            m[k] = set()
        m[k].add(v)
    for k,vs in m.items():
        m[k] = tuple(vs)
    return m

def impute(X):
    """
    Impute the value of missing genes to be the mean expression of that gene
    across all experiments.
    
    Also drops experiments and genes with no data.
    """
    if not pandas.isnull(X).sum().sum():
        return X
    X = X.dropna(axis=0, how="all").dropna(axis=1, how="all")
    model = Imputer(axis=1, strategy="mean")
    return DataFrame(model.fit_transform(X),
                     index=X.index,
                     columns=X.columns)

URSA_MEAN = 7.32
URSA_STDEV = 2.63

def symbol_map():
    m = {}
    c.execute("""
    SELECT id, symbol 
    FROM gene 
    WHERE taxon_id=9606 AND symbol IS NOT NULL""")
    return dict(c)

def collapse_to_symbols(X, axis=1, min_pct=0.2):
    # FIXME: collapse by max mean instead of max
    thresh = int(X.shape[0] * min_pct)
    return X.groupby(symbol_map(), axis=1).max().dropna(thresh=thresh, axis=1)

def fetch_genes(taxon_id):
    c.execute("""
    SELECT id, symbol, name 
    FROM gene 
    WHERE taxon_id=%s 
    ORDER BY id""", (taxon_id,))
    return DataFrame.from_items([(row[0], row) for row in c], 
                                columns=["id", "symbol", "name"],
                                orient="index")

def taxon_expression(taxon_id, require_age=True, require_gender=False, limit=None):
    c = db.cursor("taxon_expression")
    c.itersize = 1
    genes = fetch_genes(taxon_id)
    query = """
    SELECT sample.id, sample.age, sample.gender, expression.data 
    FROM expression 
    INNER JOIN sample 
    ON expression.sample_id=sample.id 
    INNER JOIN platform
    ON sample.platform_id=platform.id
    WHERE platform.taxon_id=%s"""
    if require_age:
        query += "\nAND sample.age IS NOT NULL"
    if require_gender:
        query += "\nAND sample.gender IS NOT NULL"
    if limit:
        query += "\tLIMIT " + str(limit)
    c.execute(query, (taxon_id,))
    samples, age, gender, expression = zip(*c)
    X = DataFrame.from_records(list(expression), 
                               index=samples, columns=genes.index)
    X.index.name = "Sample ID"
    X.columns.name = "Gene ID"
    P = DataFrame({"age": age, "gender": gender}, index=samples)
    P.index.name = "Sample"
    return X, P

def sample_expression(accession):
    c.execute("""
    SELECT data FROM expression
    INNER JOIN sample
    ON sample.id=expression.sample_id
    WHERE sample.accession=%s""", (accession,))
    try:
        return Series(c.__next__()[0], 
                      index=fetch_genes(9606).index)
    except StopIteration:
        return None

def infer_tissue(X):
    """
    Infer tissue for sample(s). The columns should be gene symbols.
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
    c.execute("""
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


def coo_to_df(triplets):
    """
    Create a SparseDataFrame from a sequence of (row,col,value) triplets.
    """
    data = defaultdict(dict)
    for row, col, val in triplets:
        data[col][row] = val
    return SparseDataFrame(data)

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

def ontology_graph():
    import networkx as nx
    c.execute("""
    SELECT agent, target 
    FROM term_relation
    WHERE probability=1""")
    g = nx.DiGraph()
    g.add_edges_from(c)
    return g


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

#C = X.corrwith(P)
#C = C.ix[numpy.invert(numpy.isnan(C))]
#C.sort()
#db.close()

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

def gsea_enrichment_score(labels, ranking, p=1):
    hits = np.power(np.abs(ranking * labels), p)
    p_hit = hits.cumsum() / hits.sum()
    p_miss = (np.abs(labels - 1) / (labels.shape[0] - labels.sum())).cumsum()
    dx = p_hit - p_miss
    return dx[np.abs(dx).argmax()]

import networkx as nx
from scipy.stats import fisher_exact

class Ontology(nx.DiGraph):
    def __init__(self, taxon_id=9606, namespace="GO"):
        super(Ontology, self).__init__()
        c.execute("""
        SELECT term.id, term.accession, term.name
        FROM term
        INNER JOIN ontology
        ON term.ontology_id=ontology.id
        WHERE ontology.namespace=%s
        """, (namespace,))
        for id,accession,name in c:
            self.add_node(id, accession=accession, name=name)
        c.execute("""
        SELECT agent,target
        FROM term_relation
        INNER JOIN term t1
        ON t1.id=agent
        INNER JOIN term t2
        ON t2.id=target
        INNER JOIN ontology o1
        ON t1.ontology_id=o1.id
        INNER JOIN ontology o2
        ON t2.ontology_id=o2.id
        WHERE probability=1
        AND o1.namespace=%s
        AND o2.namespace=%s""", (namespace, namespace))
        self.add_edges_from(c)
        c.execute("""
        SELECT term_id,gene_id FROM gene_term
        INNER JOIN term
        ON gene_term.term_id=term.id
        INNER JOIN ontology
        ON term.ontology_id=ontology.id
        INNER JOIN gene
        ON gene_term.gene_id=gene.id
        WHERE ontology.namespace=%s
        AND gene.taxon_id=%s
        """, (namespace, taxon_id))

        self._genes = set()

        for term_id,gene_id in c:
            self.add_node(term_id)
            data = self.node[term_id]
            if not "genes" in data:
                data["genes"] = set()
            data["genes"].add(gene_id)
            self._genes.add(gene_id)

    def enrichment(self, selected, background=None, cutoff=5e-2):
        if background is None:
            background = self._genes
        else:
            background = set(background) & self._genes

        selected = set(selected) & background

        index = []
        rows = []
        for term in self.nodes():
            annotated = self.node[term].get("genes", set())
            intersect = len(selected & annotated)
            if intersect <= 1:
                continue
            b = len(selected - annotated)
            c = len(annotated - selected)
            d = len(background) - (intersect + b + c)
            ct = [[intersect,b],[c,d]]
            odds_ratio, p_value = fisher_exact(ct, alternative="greater")
            if p_value <= cutoff:
                index.append(term)
                rows.append((term, self.node[term]["accession"],
                             self.node[term]["name"], 
                             intersect, len(annotated),
                             odds_ratio, p_value))
        columns = ["Term ID", "Accession", "Description", "Count", 
                   "Total", "Odds Ratio", "P-Value"]
        return pandas.DataFrame.from_records(rows, index="Term ID",
                                             columns=columns).sort("P-Value")

    def gsea(self, ranking, p=1, min_count=10, n_permutations=1000):
        ranking = ranking.copy()
        ranking.sort(ascending=False)
        ranked_genes = ranking.index
        ranking = numpy.array(ranking)
        rows = []

        for term in self.nodes():
            print(term)
            annotated = self.node[term].get("genes", set())
            if len(annotated) < min_count:
                continue
            labels = numpy.array([1 if g in annotated else 0 \
                                  for g in ranked_genes])
            # NOTE: actual GSEA permutes the phenotype labels on the
            # original gene expression matrix, then recomputes the
            # correlations. Since that would be really slow, I'm just
            # permuting the genes here.
            es = gsea_enrichment_score(labels, ranking)
            es_null = []
            for i in range(n_permutations):
                np.random.shuffle(labels)
                es_null.append(gsea_enrichment_score(labels, ranking))
            es_null = np.array(es_null)
            es_neg = False
            if es < 0:
                es_neg = True
                es *= -1
                es_null *= -1
            es_null = es_null[es_null > 0]
            es_null.sort()
            ix = np.searchsorted(es_null, es)
            n = es_null.shape[0]
            p_value = (n - ix) / n
            if es_neg:
                es *= -1
            rows.append((term, 
                         self.node[term]["accession"], 
                         self.node[term]["name"], 
                         labels.sum(), es, p_value))
        columns = ["Term ID", "Accession", "Description", 
                   "Count", "Enrichment Score", "P-Value"]
        return DataFrame.from_records(rows, 
                                      columns=columns, 
                                      index="Term ID").sort("P-Value")


def gene_information(taxon_id):
    c.execute("""
    SELECT id,symbol,name
    FROM gene
    WHERE taxon_id=%s""", (taxon_id,))
    return DataFrame.from_records(list(map(tuple,c)), 
                                  index="Gene ID", 
                                  columns=["Gene ID", "Symbol", "Name"])

go = Ontology()
df = go.gsea(C["Age Correlation"])

"""
X, P = taxon_expression(9606)
X = X.dropna(axis=1, how="all")
Xi = DataFrame(Imputer().fit_transform(X.as_matrix()),
               index=X.index, columns=X.columns)
"""

#C = Xi.corrwith(P["age"])
#C.name = "Age Correlation"
#C = genes.join(C, how="inner").sort("Age Correlation")

"""
ix = (P["age"] > 10) & (P["age"] < 100)
C = Xi.ix[ix,:].corrwith(P.ix[ix,"age"])
C.sort()
"""

print(C.ix[[3479,3480,2688,2690],:])
#genes = gene_information(9606)

#n = 200
#young = go.enrichment(C[:n].index, C.index)
#old = go.enrichment(C[-n:].index, C.index)
