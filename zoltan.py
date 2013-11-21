import numpy
import matplotlib
import matplotlib.pyplot as plt

import dcaf.dataset
import dcaf.learn
from dcaf.dataset import load_dataset
from dcaf.learn import QuantileNormalizer, Scaler
from dcaf.expression import aov

X, P = load_dataset("hfd_hippocampus")

def normalize_by_housekeeping_genes(X):
    # Normalize to Ywhaz/B2m/Hprt
    baseline = X.T.head(3).mean()
    return ((2 ** baseline) / (2 ** X.T)).T

X = QuantileNormalizer().fit_transform(2 ** X).apply(numpy.log2)
hk = ["Ywhaz", "B2m", "Hprt"]

#groups = list(map(lambda x,y: x+"."+y, P["Age"], P["Diet"]))
#m = aov(X, P)
#print(m.ix[hk,:])

def volcano(X, factor):
    assert len(factor) == X.shape[0]

    f1,f2 = list(set(factor))
    log10P = []
    log2FC = []

    for tx in X.columns:
        g1 = X.ix[factor==f1,tx]
        g2 = X.ix[factor==f2,tx]
        t, p = ttest_ind(g1,g2)
        fc = g1.mean() - g2.mean()
        log10P.append(- math.log10(p))
        log2FC.append(fc)

    #FIXME: may be backwards
    matplotlib.rcParams.update({"font.size": 10})
    plt.title("%s vs %s" % (f1,f2))
    plt.xlabel(r"log$_{2}$ Fold Change")
    plt.ylabel("-log$_{10}$ $p$ Value")
    plt.xlim((-3,3))
    plt.scatter(log2FC, log10P, s=40)
    for label, y, x in zip(X.columns, log10P, log2FC):
        if (abs(x) > 1) or (y > 1.75):
            plt.annotate(label, 
                         xy=(x,y), xytext=(10,10),
                         textcoords="offset points", ha="right")

#volcano(X, P["Diet"])
#volcano(X, P["Age"])
#plt.show()

# TODO: volcano plot
# TODO: Trend deviation

import dcaf.ontology
import dcaf.db

session = dcaf.db.get_session()
o = dcaf.ontology.Ontology(session)
print(o)
