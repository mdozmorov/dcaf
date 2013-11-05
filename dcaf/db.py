"""
API for SQL database containing various genomic data.
"""
#import enum
import io
import tarfile

import psycopg2

import dcaf
import dcaf.ontology
from dcaf.util import Proxy, partition
from dcaf.io import Handle

#TranscriptID = enum.Enum("Entrez", "Symbol")

class Table(object):
    def __init__(self, db, name):
        self._db = db
        self._name = name
    
    def __iter__(self):
        return self._db("SELECT * FROM %s" % self.name)
        
    @property
    def name(self):
        return self._name
        
    def columns(self):
        raise NotImplementedError
    
    def copy(self, handle, **kwargs):
        """
        Copy rows from a file handle into this table 
        using the PostgreSQL COPY command.
        """
        c = self._db.cursor()
        c.copy_from(handle, self.name, **kwargs)
        self._db.commit()
        """
        for i,group in enumerate(partition(100000, handle)):
            print(i)

            buffer = io.BytesIO()
            for line in group:
                buffer.write(line)
            buffer.seek(0)

            c.copy_from(buffer, self.name, **kwargs)
            self._db.commit()
        """
        
    def truncate(self):
        """
        Delete all the rows in this table and dependent tables.
        """
        self._db("TRUNCATE %s CASCADE" % self.name)
    
class Connection(Proxy):
    """
    Wrapper for a connection to the dcaf PostgreSQL database.
    """
    def __init__(self, user=None, password=None, 
                 database=None, host=None, port=None):
        """
        Create a new connection to the dcaf database.
        
        If the connection is successful, but the schema 
        has not yet been loaded, then the DB schema will be loaded,
        and some small, core datasets will be imported.
        
        Large datasets like MEDLINE, GEO, etc, must be imported manually.
        """
        self._db = psycopg2.connect(user=user, password=password, 
                                    database=database, 
                                    host=host, port=port)
        super(Connection,self).__init__(self._db)

        if not "taxon" in self.tables:
            self.initialize()

    def initialize(self):
        """
        Load the database schema and core datasets.
        """
        self.execute_file(dcaf.util.open_data("sql/schema.sql"))
        self.execute_file(dcaf.util.open_data("sql/functions.sql"))
        
    @property
    def tables(self):
        q = """
        SELECT table_name 
        FROM information_schema.tables
        WHERE table_schema = 'public'
        """
        return list(row[0] for row in self.execute(q))
        
    def table(self, name):
        return Table(self, name)
        
    def execute(self, sql, *params):
        """
        Execute SQL on this connection.
        
        :param sql: Raw SQL to be executed
        :type sql: str
        :param params: Parameters to be replaced in the format string
        :returns: A cursor object containing any applicable result set.
        """
        c = self._db.cursor()
        c.execute(sql, params)
        self._db.commit()
        return c

    def __call__(self, sql, *params):
        """
        Shorthand for Connection.execute.
        """
        return self.execute(sql, *params)

    def execute_file(self, file):
        """
        Execute SQL contained in a file handle or file path.

        :param file: File path or file-like object
        :type file: str or handle
        :returns: A cursor object containing any applicable result set.
        """
        if not hasattr(file, "read"):
            file = open(file)
        with file:
            return self.execute(file.read())
    
    def ontology(self, taxon_id=9606, namespace="GO"):
        """
        Load an ontology from this database into a directed graph.
        """
        return dcaf.ontology.Ontology(self._db.cursor(), 
                                      taxon_id=taxon_id, namespace=namespace)
        
    def genes(self, taxon_id):
        c = self("""
        SELECT id, symbol, name 
        FROM gene 
        WHERE taxon_id=%s 
        ORDER BY id""", (taxon_id,))
        return DataFrame.from_items([(row[0], row) for row in c], 
                                    columns=["id", "symbol", "name"],
                                    orient="index")

    def expression(self, taxon_id=9606, 
                   id_type="entrez",
                   require_age=True, 
                   require_gender=False, 
                   limit=None):
        """
        Get an expression set corresponding to the given criteria.
        """
        assert(id_type in ("entrez", "symbol"))

        genes = self.genes(taxon_id)
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

        c = self._db.cursor("expression")
        c.itersize = 1
        c.execute(query, (taxon_id,))
        samples, age, gender, expression = zip(*c)
        X = DataFrame.from_records(list(expression), 
                                   index=samples, columns=genes.index)
        X.index.name = "Sample ID"
        X.columns.name = "Gene ID"
        P = DataFrame({"age": age, "gender": gender}, index=samples)
        P.index.name = "Sample"

        if id_type == "symbol":
            X = self._collapse_to_symbols(X, taxon_id)

        return X, P

    def sample_expression(self, accession):
        """
        Get expression for a single sample.
        """
        c = self("""
        SELECT data FROM expression
        INNER JOIN sample
        ON sample.id=expression.sample_id
        WHERE sample.accession=%s""", (accession,))
        try:
            return Series(c.__next__()[0], 
                          index=self.genes(9606).index)
        except StopIteration:
            return None

    def _symbol_map(self, taxon_id):
        m = {}
        return dict(self("""
        SELECT id, symbol 
        FROM gene 
        WHERE taxon_id=%s AND symbol IS NOT NULL""", taxon_id))

    def _collapse_to_symbols(self, X, taxon_id, axis=1, min_pct=0.2):
        # FIXME: collapse by max mean instead of max
        thresh = int(X.shape[0] * min_pct)
        return X.groupby(self._symbol_map(taxon_id), 
                         axis=1).max().dropna(thresh=thresh, axis=1)


if __name__ == "__main__":
    import sys
    BASE = "ftp://ftp.ncbi.nlm.nih.gov"

    db = Connection(database="testdb")

    #print("Loading taxon ...")
    #db.table("taxon").copy(
    #    Handle(BASE + "/pub/taxonomy/taxdump.tar.gz") \
    #    | "tar Oxf - names.dmp | grep scientific | cut -f1,3")

    print("Loading gene ...")
    db.table("gene").truncate()
    db.table("gene").copy(
    Handle(BASE + "/gene/DATA/gene_info.gz") \
        | "sed 1d | grep -v '^#' | cut -f1-3,12 | awk '$1 != 366646'",
        null="-")

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

def gsea_enrichment_score(labels, ranking, p=1):
    hits = np.power(np.abs(ranking * labels), p)
    p_hit = hits.cumsum() / hits.sum()
    p_miss = (np.abs(labels - 1) / (labels.shape[0] - labels.sum())).cumsum()
    dx = p_hit - p_miss
    return dx[np.abs(dx).argmax()]

def gene_information(taxon_id):
    c.execute("""
    SELECT id,symbol,name
    FROM gene
    WHERE taxon_id=%s""", (taxon_id,))
    return DataFrame.from_records(list(map(tuple,c)), 
                                  index="Gene ID", 
                                  columns=["Gene ID", "Symbol", "Name"])

#go = Ontology()
#df = go.gsea(C["Age Correlation"])

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

#print(C.ix[[3479,3480,2688,2690],:])
#genes = gene_information(9606)

#n = 200
#young = go.enrichment(C[:n].index, C.index)

def queryUCSC(assembly, query):
    db = mysql.connector.connect(user="genome", password="genome",
                                 host="wrendb", database=assembly)
    df = pandas.io.sql.read_frame(query, db)
    db.close()
    return df
