import networkx as nx
import pandas
import numpy

from scipy.stats import fisher_exact
from pandas import DataFrame

from .db.model import *
from .db.model import Ontology as DBOntology

def gsea_enrichment_score(labels, ranking, p=1):
    hits = numpy.power(numpy.abs(ranking * labels), p)
    p_hit = hits.cumsum() / hits.sum()
    p_miss = (numpy.abs(labels - 1) / (labels.shape[0] - labels.sum())).cumsum()
    dx = p_hit - p_miss
    return dx[numpy.abs(dx).argmax()]

class Ontology(nx.DiGraph):
    def __init__(self, session, taxon_id=9606, namespace="GO"):
        super(Ontology, self).__init__()
        o = session.query(DBOntology).filter_by(namespace=namespace)

        for t in session.query(Term).join(DBOntology.terms).filter(DBOntology.namespace==namespace):
            self.add_node(t.id, accession=t.accession, name=t.name)

        # FIXME: figure out multiple joins in SQLAlchemy
        for e in session.query(TermRelation):
            if e.agent in self and e.target in self and e.probability == 1:
                self.add_edge(e.agent, e.target)

        self._genes = set()

        for e in session.query(GeneTerm).join(Term).join(DBOntology).join(Gene)\
            .filter(DBOntology.namespace==namespace).filter(Gene.taxon_id==taxon_id):
            self.add_node(e.term_id)
            data = self.node[e.term_id]
            if not "genes" in data:
                data["genes"] = set()
            data["genes"].add(e.gene_id)
            self._genes.add(e.gene_id)

    def enrichment(self, selected, background=None, cutoff=5e-2):
        if background is None:
            background = self._genes
        else:
            background = set(background) & self._genes

        selected = set(selected) & background

        index = []
        rows = []
        for term in self.nodes():
            annotated = self.node[term].get("genes", set()) & background
            intersect = len(selected & annotated)
            if intersect <= 1:
                continue
            b = len(selected - annotated)
            c = len(annotated - selected)
            d = len(background) - (intersect + b + c)
            ct = [[intersect,b],[c,d]]
            assert (intersect + b + c + d) == len(background)
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
                numpy.random.shuffle(labels)
                es_null.append(gsea_enrichment_score(labels, ranking))
            es_null = numpy.array(es_null)
            es_neg = False
            if es < 0:
                es_neg = True
                es *= -1
                es_null *= -1
            es_null = es_null[es_null > 0]
            es_null.sort()
            ix = numpy.searchsorted(es_null, es)
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

if __name__ == "__main__":
    #go = Ontology()
    #df = go.gsea(C["Age Correlation"])

    #print(C.ix[[3479,3480,2688,2690],:])
    #genes = gene_information(9606)

    #n = 200
    #young = go.enrichment(C[:n].index, C.index)
    pass
