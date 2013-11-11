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
                if len(self._handles) == 3:
                    break
        
        self._global_mean = numpy.array([h.summary.sumData / h.summary.basesCovered
                                         for h in self._handles])
        self._index = [os.path.basename(h._path).rstrip(".bw") 
                       for h in self._handles]

    def count_matrix(self, regions, normalize=True):
        rows = {}
        for r in regions:
            rows[r.name] = self.count(r, normalize=normalize)
        return pandas.DataFrame(rows)
        
    def count(self, region, normalize=False):
        r = region
        data = numpy.array([h.summarize_region(r.chrom, r.start, r.end).mean
                            for h in self._handles])
        if normalize:
            data /= self._global_mean
        return pandas.Series(data, index=self._index)
        
if __name__ == "__main__":
    from itertools import islice

    bws = BigWigSet("~/data/continuous-tracks/hg19/RNAseq/")

    with BEDFile("~/data/knownGene.bed.gz") as h:
        regions = list(islice(h, 10))
    df = bws.count_matrix(regions)
    print(df.corrwith(df.iloc[:,0]))
    print(df)
