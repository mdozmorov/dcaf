"""
Various types of indices for genomic intervals, both RAM and disk-based.
"""

from dcaf.data_structure import IntervalTree

__all__ = ["GenomicRegionRAMIndex"]

class GenomicRegionRAMIndex(object):
    """
    A RAM-based index for genomic regions on multiple chromosomes.
    
    Internally, it stores one IntervalTree per chromosome.
    """
    def __init__(self):
        self._intervals = {}
        self._roots = {}
        self._built = False

    def add(self, chrom, start, end, data=None):
        """
        Add another interval to the index. This method can only be
        called if the index has not been built yet.
        """
        if self._built:
            raise Exception("Can't currently mutate a constructed IntervalTree.")
        self._intervals.setdefault(chrom, [])
        self._intervals[chrom].append((start, end, data))
    
    def build(self):
        """
        Build the index. After this method is called, new intervals
        cannot be added.
        """
        assert(self._intervals)
        for chrom, intervals in self._intervals.items():
            self._roots[chrom] = IntervalTree(intervals)
        self._built = True
    
    def search(self, chrom, start, end):
        """
        Search the index for intervals overlapping the given 
        chrom, start, and end.
        """
        if not self._built:
            raise Exception("Must call %s.build() before using search()" % \
                            str(type(self)))
        return self._roots[chrom].search(start, end)
    
    def __iter__(self):
        assert(self._built)
        for chrom, node in self._roots.items():
            yield from node
