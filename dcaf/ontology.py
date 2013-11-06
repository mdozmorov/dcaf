import networkx as nx
from scipy.stats import fisher_exact

class Ontology(nx.DiGraph):
    def __init__(self, cursor, taxon_id=9606, namespace="GO"):
        super(Ontology, self).__init__()
        c = cursor
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

def gsea_enrichment_score(labels, ranking, p=1):
    hits = np.power(np.abs(ranking * labels), p)
    p_hit = hits.cumsum() / hits.sum()
    p_miss = (np.abs(labels - 1) / (labels.shape[0] - labels.sum())).cumsum()
    dx = p_hit - p_miss
    return dx[np.abs(dx).argmax()]

#go = Ontology()
#df = go.gsea(C["Age Correlation"])

#print(C.ix[[3479,3480,2688,2690],:])
#genes = gene_information(9606)

#n = 200
#young = go.enrichment(C[:n].index, C.index)
