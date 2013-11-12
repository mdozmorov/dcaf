"""
Stuff that's still in an early stage of development. Use at your own risk.
"""
import os
import gzip

from collections import namedtuple

import numpy
import pandas

BED = namedtuple("BED", "chrom start end name score strand rest")

class BEDFile(object):
    # FIXME: ignore genome browser fields

    def __init__(self, path):
        path = os.path.expanduser(path)
        self._handle = gzip.open(path, "rt")

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()
    
    def close(self):
        self._handle.close()
    
    def __iter__(self):
        for line in self._handle:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            try:
                chrom = fields[0]
                start = int(fields[1])
                end = int(fields[2])
                name = fields[3] if len(fields) > 3 else None
                score = None
                if len(fields) > 4:
                    score = float(fields[4]) if fields[4] != "." else numpy.nan
                strand = fields[5] if len(fields) > 5 else None
                yield BED(chrom, start, end, name, score, strand, fields[6:])
            except Exception as e:
                raise ValueError("Malformed BED line:\n" + line)

from dcaf.io import BigWigFile

def _count_regions(args):
    path, regions = args
    h = BigWigFile(path)
    counts = {}
    for r in regions:
        mu = h.summarize_region(r.chrom, r.start, r.end).mean
        mu *= h.summary.basesCovered / h.summary.sumData
        counts[r.name] = mu
    return pandas.Series(counts)

class BigWigSet(object):
    def __init__(self, basedir, recursive=False):
        # FIXME: actually test the file magic number to see if it is a BigWig
        # FIXME: implement recursive search
        # FIXME: sort paths and return __hash__ as hash of paths
        # FIXME (?): incorporate standard deviation into normalization
        # TODO: implement lazy/memoized index loading in BBI files
        # TODO: make count_matrix parallel
        assert not recursive

        self._handles = []
        basedir = os.path.expanduser(basedir)
        for file in os.listdir(basedir):
            if file.endswith(".bw"):
                path = os.path.join(basedir, file)
                self._handles.append(BigWigFile(path))
        
        self._global_mean = numpy.array([h.summary.sumData / h.summary.basesCovered
                                         for h in self._handles])
        self._index = [os.path.basename(h._path).rstrip(".bw") 
                       for h in self._handles]

    def count_matrix(self, regions, normalize=True):
        import multiprocessing as mp
        p = mp.Pool(mp.cpu_count())
        rows = list(p.map(_count_regions, 
                          ((h._path, regions) for h in self._handles)))
        return pandas.DataFrame(rows, index=self._index)
          
    def count(self, region, normalize=False):
        r = region
        data = numpy.array([h.summarize_region(r.chrom, r.start, r.end).mean
                            for h in self._handles])
        if normalize:
            data /= self._global_mean
        return pandas.Series(data, index=self._index)
        
def entrez_gene_loci():
    import dcaf.db
    from itertools import groupby
    import operator

    db = dcaf.db.UCSCConnection.from_configuration()
    q = """
    SELECT DISTINCT chrom,txStart,txEnd,
      CAST(knownToLocusLink.value AS INTEGER) AS geneID,'.',strand
    FROM knownGene
    INNER JOIN knownToLocusLink
    ON knownGene.name=knownToLocusLink.name
    ORDER BY chrom,geneID,txStart,txEnd
    """
    for name, rows in groupby(db(q), operator.itemgetter(3)):
        prev = BED(*next(rows) + (None,))
        for item in rows:
            item = BED(*item + (None,))
            if item.start <= prev.end:
                prev = BED(prev.chrom, prev.start, item.end, 
                           name, prev.score, prev.strand, None)
            else:
                yield prev
                prev = item
        yield prev

def main():
    import sys
    from itertools import islice

    #regions = list(islice(entrez_gene_loci(), 10))
    regions = list(entrez_gene_loci())
    #print(_count_regions("/home/gilesc/data/RNAseq/bw/DRR001622.bw", regions))

    #with BEDFile("~/data/knownGene.bed.gz") as h:
    #    regions = list(islice(h, 10))

    bws = BigWigSet("~/data/RNAseq/bw/")
    df = bws.count_matrix(regions)
    #import sys
    df.to_csv(sys.stdout, sep="\t", float_format="%0.3f")

    #print(df.corrwith(df.iloc[:,0]))
    #print(df)

if __name__ == "__main__":
    import dcaf.io
    bws = BigWigSet("~/data/RNAseq/bw/")
    r = BED("chr7", 13141016, 13743774, None, None, None, None)
    v = bws.count(r, normalize=True)
    X = dcaf.io.read_matrix("hg19.eg.mat")
    c = X.corrwith(v)
    c.sort(ascending=False)
    c.to_csv("AC011288.2.cor", sep="\t")
