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
