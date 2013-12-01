from itertools import repeat

import numpy
import pandas
import matplotlib
import matplotlib.pyplot as plt

from statsmodels.stats.multitest import multipletests
from scipy.stats import ttest_ind
from patsy import dmatrix
from statsmodels.regression.linear_model import OLS

import dcaf.dataset
import dcaf.learn
from dcaf.dataset import load_dataset
from dcaf.learn import QuantileNormalizer, Scaler
from dcaf.expression import aov
from dcaf.db.model import Gene

def normalize_by_housekeeping_genes(X):
    # Normalize to Ywhaz/B2m/Hprt
    baseline = X.T.head(3).mean()
    return ((2 ** baseline) / (2 ** X.T)).T

def overlaps(s1,e1,s2,e2):
    return (s2 < e1) and (s1 < e2)

def remove_point_overlaps(points, shape=(0.5,0.25)):
    # Orders rectangles by distance from (0,0), then
    # removes any that overlap with more distant 
    # rectangles

    # FIXME: implement RTree
    points.sort(key=lambda pt: -(abs(pt[0]) + abs(pt[1])))
    dx,dy = numpy.array(shape) / 2
    o_points = []
    for x1,y1,label in points:
        for x2,y2,_ in o_points:
            x1_s, x1_e = (x1-dx),(x1+dx)
            x2_s, x2_e = (x2-dx),(x2+dx)
            y1_s, y1_e = (y1-dy),(y1+dy)
            y2_s, y2_e = (y2-dy),(y2+dy)
            if overlaps(x1_s,x1_e,x2_s,x2_e) and \
               overlaps(y1_s,y1_e,y2_s,y2_e):
                break
        else:
            o_points.append((x1,y1,label))
    return o_points

"""
def volcano1(X, factor):
    assert len(factor) == X.shape[0]

    f1,f2 = list(set(factor))
    log10P = []
    log2FC = []

    for tx in X.columns:
        g1 = 2 ** X.ix[factor==f1,tx]
        g2 = 2 ** X.ix[factor==f2,tx]
        t, p = ttest_ind(g1,g2)
        fc = math.log2(g1.mean()) - math.log2(g2.mean())
        log10P.append(- math.log10(p))
        log2FC.append(fc)

    #FIXME: may be backwards
"""

def volcano(fold_change, p_value, labels,
            title=None):
    """
    Make a volcano plot for visualizing fold changes
    and p-values.
    """

    matplotlib.rcParams.update(
        {
            "font.size": 12,
            "xtick.labelsize": 15,
            "ytick.labelsize": 15,
            "xtick.major.size": 10,
            "ytick.major.size": 10
        })

    if title:
        plt.title(title, size=28)

    plt.xlabel(r"log$_{2}$ Fold Change", size=20)
    plt.ylabel(r"-log$_{10}$ $p$ Value", size=20)
    x_max = numpy.ceil(numpy.fabs(fold_change)).max()
    plt.xlim((-x_max,x_max))
    plt.scatter(fold_change, p_value, s=40)
    
    points = []
    for label, y, x in zip(labels, p_value, fold_change):
        if (abs(x) > 0.58) or (y > 1.75):
            points.append((x,y,label))
    
    points = remove_point_overlaps(points)
    
    for x,y,label in points:
        plt.annotate(label, 
                     xy=(x,y), xytext=(10,10),
                     textcoords="offset points", 
                     ha="right")

def aov_plot_volcano(aov, factor):
    volcano(aov[factor]["coef"],
            -numpy.log10(aov[factor]["p"]),
            aov.index,
            title=factor)

def multiaov(X,D,correct=True):
    data = {}
    for tx in X.columns:
        model = OLS(X.ix[:,tx], exog=D).fit()
        data[tx] = list(model.params) + \
                   list(model.pvalues)
    factors = model.params.index
    columns = list(zip(factors, repeat("coef"))) + \
              list(zip(factors, repeat("p")))
    aov = pandas.DataFrame.from_dict(data).T
    aov.index.names = ("Transcript",)
    aov.columns = pandas.MultiIndex.from_tuples(columns)
    aov.columns.names = ("Factor", "Metric")
    
    if correct:
        for ix in aov.columns:
            if ix[1] != "p":
                continue
            aov[ix] = multipletests(aov[ix],
                                    method="h")[1]
    return aov


X, P = load_dataset("hfd_hippocampus")
X = QuantileNormalizer().fit_transform(2 ** X)\
                        .apply(numpy.log2)

f_age = "C(Age,levels=['Young','Old'])"
f_diet = "C(Diet,levels=['Control','HFD'])"
f = " + ".join([f_age, f_diet, f_age+":"+f_diet])
D = dmatrix(f, data=P, return_type="dataframe")
D.columns = ["Intercept", "Age", "Diet", "Age:Diet"]

X_g = X.groupby([P["Age"],P["Diet"]]).mean()

#aov = multiaov(X,D,correct=False)
#aov_plot_volcano(aov, "Age:Diet")
#plt.show()

# TODO: Trend deviation

import imp

import dcaf.ontology
imp.reload(dcaf.ontology)
import dcaf.db

session = dcaf.db.get_session()

symbols = {}
for gene in session.query(Gene).filter_by(taxon_id=10090):
    symbols[gene.symbol] = gene.id

o = dcaf.ontology.Ontology(session, taxon_id=10090)

X_eg = X.T.groupby(symbols).max().T

aov = multiaov(X_eg, D)
#p = list(aov[("Age:Diet", "p")].order().index)
#ea = o.enrichment(p[:10], background=p, cutoff=1)
gsea = o.gsea(aov[("Age:Diet", "p")])
