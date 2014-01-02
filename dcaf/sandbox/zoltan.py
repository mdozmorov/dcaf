import os
import imp
from itertools import repeat

import numpy
import pandas
import matplotlib
import matplotlib.pyplot as plt
from statsmodels.stats.multitest import multipletests
from scipy.stats import ttest_ind
from patsy import dmatrix
from statsmodels.regression.linear_model import OLS
from openpyxl import Workbook

import dcaf.dataset
import dcaf.learn
import dcaf.db
import dcaf.expression
from dcaf.ontology import Ontology
from dcaf.dataset import load_dataset
from dcaf.learn import QuantileNormalizer, Scaler
from dcaf.db.model import Gene

# http://www.gene-quantification.de/livak-2001.pdf

RESET = "X_g" not in locals()
#RESET = True

OUTDIR = "/tmp/zoltan/"
if not os.path.exists(OUTDIR):
    os.makedirs(OUTDIR)

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

def volcano(fold_change, p_value, labels,
            title=None, xlabel=None, ylabel=None):
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
    xlabel = xlabel or r"log$_{2}$ Fold Change"
    ylabel = ylabel or r"-log$_{10}$ $p$ Value"

    plt.xlabel(xlabel, size=20)
    plt.ylabel(ylabel, size=20)
    x_max = numpy.ceil(numpy.fabs(fold_change) + 0.5).max()
    plt.xlim((-x_max,x_max))
    plt.scatter(fold_change, p_value, s=40)
    
    points = list(zip(fold_change, p_value, labels))
    points.sort(key=lambda pt: -(abs(pt[0])*pt[1]))
    points = points[:20]
    points = remove_point_overlaps(points)
    points = points[:10]
    
    for x,y,label in points:
        plt.annotate(label, 
                     xy=(x,y), xytext=(10,10),
                     textcoords="offset points", 
                     ha="right")

def aov_plot_volcano(aov, factor, **kwargs):
    if not "title" in kwargs:
        kwargs["title"] = factor

    volcano(aov[factor]["coef"],
            -numpy.log10(aov[factor]["p"]),
            aov.index, **kwargs)

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


if RESET:
    X, P = load_dataset("hfd_hippocampus")
    #X = QuantileNormalizer().fit_transform(2 ** X)\
    #                        .apply(numpy.log2)
    X = normalize_by_housekeeping_genes(X)

    f_age = "C(Age,levels=['Young','Old'])"
    f_diet = "C(Diet,levels=['Control','HFD'])"
    f = " + ".join([f_age, f_diet, f_age+":"+f_diet])
    D = dmatrix(f, data=P, return_type="dataframe")
    D.columns = ["Intercept", "Age", "Diet", "Age:Diet"]

    X_g = X.groupby([P["Age"],P["Diet"]]).mean()

if RESET:

    for correct in [True, False]:
        aov = multiaov(X,D,correct=correct)
        ext = ("_corrected" if correct else "_nominal") + ".png"

        FIG_KWARGS = {"pad_inches":0.5, "bbox_inches":"tight", "dpi":90}

        aov_plot_volcano(aov, "Age")
        plt.savefig(OUTDIR+"Age_volcano" + ext, **FIG_KWARGS)
        plt.clf()

        aov_plot_volcano(aov, "Diet")
        plt.savefig(OUTDIR+"Diet_volcano" + ext, **FIG_KWARGS)
        plt.clf()

        aov_plot_volcano(aov, "Age:Diet", 
                         title="Age-Diet Interaction",
                         xlabel=r"ANOVA Coefficient (log$_{2}$ Fold Change)")
        plt.savefig(OUTDIR+"AgeDietInteraction_volcano" + ext, **FIG_KWARGS)
        plt.clf()


if RESET:
    session = dcaf.db.get_session()

    symbols = {}
    for gene in session.query(Gene).filter_by(taxon_id=10090):
        symbols[gene.symbol] = gene.id

    X_eg = X.T.groupby(symbols).max().T

    ontologies = {}
    ontologies["MSigDB"] = Ontology(session, taxon_id=10090, namespace="MSigDB")
    ontologies["GO"] = Ontology(session, taxon_id=10090, namespace="GO")
    
def fix_column_widths(path):
    """
    Adjust the column widths in an Excel 2007 (.xslx) file 
    to the widest cell in each column.
    """
    import openpyxl
    wb = openpyxl.load_workbook(path)
    for sheet in wb.worksheets:
        col_index_names = list(sorted(sheet.column_dimensions.keys(),
                                      key = lambda k: (len(k), k)))
        for j in range(sheet.get_highest_column()):
            col = sheet.columns[j]
            col_index_name = col_index_names[j]
            width = min(max(len(str(cell.value)) for cell in col) + 0.1, 50)
            sheet.column_dimensions[col_index_name].width = width
    wb.save(path)

if RESET:
    groups = ["Age", "Diet", "Age:Diet"]
    aov = multiaov(X_eg, D, correct=False)
    gsea = []

    for group in groups:
        for o_name, o in ontologies.items():
            p = list(aov[(group, "p")].order().index)
            rs = o.gsea(aov[(group, "p")], min_count=3)
            gsea.append((group, o_name, rs))

    symbols = dict((g.id, g.symbol) for g in 
                   session.query(Gene).filter_by(taxon_id=10090))

    from pandas import ExcelWriter
    path = OUTDIR+"hfd_hippocampus.xlsx"
    writer = ExcelWriter(path)

    aov = multiaov(X_eg, D, correct=False)
    aov.index.name = "Entrez ID"
    aov[("Symbol", None)] = [symbols.get(entrez_id) for entrez_id in 
                             aov.index]
    aov = aov[[("Symbol", None)] + list(aov.columns[:-1])]
    aov.to_excel(writer, sheet_name="DE Genes", float_format="%0.3f")

    for group, name, df in gsea:
        sheet_name = (group + "_" + name).replace(":", "-")
        df.to_excel(writer,sheet_name=sheet_name, 
                    index=False, float_format="%0.3f")
    writer.save()
    fix_column_widths(path)

print("Done")
# TODO: Trend deviation
