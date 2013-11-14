import numpy

from sklearn.datasets import make_classification
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import Imputer
from sklearn.metrics import f1_score
from sklearn.svm import SVC

from dcaf.learn import InvalidColumnRemover, QuantileNormalizer

def add_random_nan(X, percent=0.3):
    ix = numpy.arange(X.size)
    numpy.random.shuffle(ix)
    X = X.copy()
    X.flat[ix[:int(percent * X.size)]] = numpy.nan
    return X

def test_invalid_column_remover():
    X,y = make_classification(1000, 100)
    X = add_random_nan(X)

    model = InvalidColumnRemover(min_count=700)
    columns_removed = X.shape[1] - model.fit_transform(X).shape[1]
    assert (columns_removed > 0) and (columns_removed < 100)

"""
def test_pipeline_improves_performance():
    X,y = make_classification(1000, 100)
    #X *= numpy.random.normal(loc=10, scale=100, size=X.shape[1])
    #print(X)

    model = SVC()

    p = Pipeline([
        ('rm-invalid', InvalidColumnRemover()),
        ('impute', Imputer()),
        #('qnorm', QuantileNormalizer()),
        ('svm', SVC())
    ])

    
    X = add_random_nan(X)

    p.fit(X,y)
    pipeline_score = f1_score(p.predict(X), y)
    
    X.flat[numpy.isnan(X.flat)] = 0
    model.fit(X,y)
    bare_model_score = f1_score(model.predict(X), y)

    print()
    print(pipeline_score, bare_model_score)

    assert pipeline_score > bare_model_score
"""
