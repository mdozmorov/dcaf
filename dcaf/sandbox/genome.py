""" 
Generic manipulation of genomic intervals.

Manipulate genomic interval sets generically. Similar in concept
to pybedtools, except it operates in memory whenever possible rather
than using tempfiles and subprocesses, and it works in Python 3. 
"""

import gzip 
import io
import os 
import sys
import warnings

from itertools import groupby, chain, islice
from operator import attrgetter
from functools import reduce, lru_cache
from collections import Iterator, Sized

import pysam
import numpy

import dcaf.util

def ensure_iter(seq):
    return seq if isinstance(seq, Iterator) else iter(seq)

def read_lines(path):
    if path.endswith(".gz"):
        return (line.decode("ascii") for line in gzip.open(path))
    else:
        return io.open(path)

class Interval(object):
    __slots__ = ["contig", "start", "end", "name", "score", "strand"]

    def __init__(self, contig, start, end=None, name=None, score=None, strand=None):
        self.contig = sys.intern(contig) if contig is not None else None
        self.start = start
        self.end = start + 1 if end is None else end
        self.name = sys.intern(name) if name is not None else None
        self.score = score
        self.strand = sys.intern(strand) if strand is not None else None
            
    def __len__(self):
        return self.end - self.start

    def __str__(self):
        def fields():
            return (getattr(self, s) for s in Interval.__slots__)
        return "\t".join([str(field) if field is not None else "." for field in fields()])
    
    def __repr__(self):
        return "<Interval %s:%d-%d>" % (self.contig, self.start, self.end)
    
    def __lt__(self, o):
        return (self.contig, self.start, self.end) < (o.contig, o.start, o.end)
    
    def __eq__(self, o):
        for slot in Interval.__slots__:
            this = getattr(self, slot)
            other = getattr(o, slot)
            if (this is not None) or (other is not None):
                if this != other:
                    return False
        return True

    def intersect(self, o):
        if self.distance(o) <= 0:
            return Interval(self.contig, 
                            max(self.start, o.start), min(self.end, o.end))

    def distance(self, o):
        """
        Return the distance in bases between this Interval and another Interval.
        Will be 0 if they have the same start and end, negative if they overlap,
        and positive if they don't overlap but are on the same contig.
        Returns infinity if they aren't on the same contig.
        """
        return max(self.start, o.start) - min(self.end, o.end) \
            if self.contig == o.contig else numpy.inf
    
    def overlaps(self, o):
        """
        Returns a boolean representing whether this Interval overlaps
        the given Interval.
        """
        return self.distance(o) <= 0
    
    @property
    def midpoint(self):
        return (self.end + self.start) / 2
    
    @staticmethod
    def parse_line(line, type="bed"):
        """
        Parse a line from a flat-file into an Interval.

        Available types:
        * bed
        """
        fields = line.strip().split("\t")
        return Interval(fields[0], int(fields[1]), int(fields[2]),
                        fields[3] if (len(fields) > 3 and fields[3] != ".") else None,
                        float(fields[4]) 
                        if (len(fields) > 4 and fields[4] != ".") else None,
                        fields[5] 
                        if (len(fields) > 5 and fields[5] != ".") else None)

class Intervals(Iterator, Sized):
    """
    An iterator containing Interval objects, with extra attributes
    and methods.
    """
    __slots__ = ["_intervals", "is_sorted"]

    def __init__(self, intervals, is_sorted=False):
        self._intervals = ensure_iter(intervals)
        self.is_sorted = is_sorted
    
    def __iter__(self):
        return self 
    
    def __next__(self):
        return next(self._intervals)

    def __len__(self):
        return reduce(lambda acc, item: acc + 1, self, 0)
    
    def save(self, path_or_handle):
        """
        Save this set of intervals to a path or file-like handle.
        
        If a path is provided, will return an IntervalFile referring to
        the written interval file. If a handle is provided, returns None.
        """
        if isinstance(path_or_handle, str):
            with open(path_or_handle, "w") as handle:
                self.save(handle)
            return IntervalFile(path_or_handle)
        elif hasattr(path_or_handle, "write"):
            for interval in self:
                print(str(interval), file=path_or_handle)
        else:
            raise Exception("Argument is not a string or writable file-like object.")
    
    def print(self):
        """
        Print this interval collection to stdout. 

        Convenience shortcut for Intervals.save(sys.stdout).
        """
        self.save(sys.stdout)

    def sort(self):
        """
        Sort intervals by contig, then start, then end.
        """
        if not self.is_sorted:
            self._intervals = iter(sorted(self._intervals))
            self.is_sorted = True
        return self

    def merge(self):
        """
        Merge adjacent overlapping intervals.
        """
        self.sort()
        def merge():
            current = next(self._intervals)
            for interval in self._intervals:
                if interval.overlaps(current):
                    current.start = min(current.start, interval.start)
                    current.end = max(current.end, interval.end)
                else:
                    yield current
                    current = interval
            yield current
        return Intervals(merge(), is_sorted=True)
    
    def uniq(self):
        """
        Delete any duplicate intervals in this collection.
        """
        self.sort()
        def uniq():
            current = next(self)
            for interval in self:
                if current != interval:
                    yield current
                current = interval
            yield current
        return Intervals(uniq(), is_sorted=True)

    def midpoints(self, unique=True):
        """
        Return a dict mapping contigs to a sorted numpy array of midpoints
        corresponding to the intervals in this set.
        
        If unique is True, the midpoints are deduplicated.
        """
        midpoints = {}
        for contig, intervals in groupby(self.sort(), attrgetter("contig")):
            midpoints[contig] = numpy.array(list(map(attrgetter("midpoint"), intervals)))
            if unique:
                midpoints[contig] = numpy.unique(midpoints[contig])
        return midpoints
   
    def head(self, n=10):
        def head():
            for i, interval in enumerate(self._intervals):
                if i < n:
                    yield interval
                else:
                    break
        return Intervals(head(), is_sorted=self.is_sorted)

class IntervalPairs(object):
    """
    An iterator over pairs (2-tuples) of Interval objects.
    """
    def __init__(self, pairs):
        self._pairs = ensure_iter(pairs)
    
    def left(self):
        return Intervals(pair[0] for pair in self._pairs)
    
    def right(self):
        return Intervals(pair[1] for pair in self._pairs)
        
    def consolidate(self):
        """
        For each pair in this set of Interval pairs, consolidate
        the pair to the region of intersection between the two regions. 
        """
        def consolidate():
            for left, right in self:
                intersection = left.intersect(right)
                if intersection:
                    yield intersection
        return Intervals(consolidate())
        
    def distance(self):
        """
        For each pair, return the distance between the two intervals.
        """
        for left,right in self:
            yield left.distance(right)

    def __iter__(self):
        return self

    def __next__(self):
        return next(self._pairs)

class IntervalFile(object):
    """
    A reader for various file formats containing genomic intervals.
    """
    def __init__(self, path):
        self._path = path
       
    def __iter__(self):
        return Intervals(Interval.parse_line(line, type="bed") 
                         for line in read_lines(self._path))
    
    def __len__(self):
        return len(iter(self))

    def __getattr__(self, attr):
        return getattr(iter(self), attr)
   
class Index(object):
    """
    An abstract query interface for indexed genomic interval sets.
    
    Inheritors must implement:
        Methods:
            fetch
            __iter__
    
    Inheritors may override:
        Fields:
            _genome (a dict of contig names to lengths)
    """

    def fetch(self, query):
        raise NotImplementedError
    
    def __iter__(self):
        raise NotImplementedError

    @property
    def _genome(self):
        """
        In the case where the genome (contig lengths) is not known, 
        it is approximated by iterating through all the intervals
        and taking the largest endpoint on each contig as the contig length.
        
        Since this is a time-consuming operation, the result is memoized.
        """
        if not hasattr(self, "__genome"):
            self.__genome = {}
            for contig, intervals in groupby(self, attrgetter("contig")):
                self.__genome[contig] = max(i.end for i in intervals)
        return self.__genome

    def intersect(self, iset):
        """
        Return a generator of interval pairs which consist of the overlapping
        regions between this index and the query interval set.
        """
        def intersect():
            for q in iset:
                for i in self.fetch(q):
                    yield (i,q)
        return IntervalPairs(intersect())

    def closest(self, item, n=1, direction="both"):
        """
        Accepts either an Interval or an Intervals as argument.
        
        If an Interval is provided:
            Return an iterable containing all the intervals from this Index on the same
            contig as the query interval, sorted by increasing distance 
            from the query interval. 
        
        If an Intervals is provided:
            Return an IntervalPairs such that, for each interval in the given Intervals, 
            the n closest intervals in this index are returned paired with the query interval.
        """
        raise NotImplementedError
        if isinstance(item, Interval):
            return Intervals(self._closest(item, direction))
        else:
            return IntervalPairs((i,c) for i in item for c in self.closest(i,n=n,direction=direction).head(n))

    def _closest(self, query, direction="both", key="midpoint"):
        contig = query.contig

        if key == "midpoint":
            keyfn = lambda x: abs(query.midpoint - x.midpoint)
            start = end = query.midpoint
        elif key == "distance":
            keyfn = query.distance
            start = query.start
            end = query.end
            for i in sorted(self.fetch(query), key=keyfn):
                yield i
        else:
            raise Exception

        step = 1000
        contig_len = self._genome[contig]
 
        while (start > 0) or (end < contig_len):
            intervals = []
            start = max(0, start - step)
            end = min(contig_len, end + step)
            if direction in ("both", "less") and (start > 0):
                q = Interval(contig, start, start + step)
                intervals.extend(self.fetch(q))
            if direction in ("both", "greater") and (end < contig_len):
                q = Interval(contig, end - step, end)
                intervals.extend(self.fetch(q))
            if not intervals:
                step *= 2
            for i in sorted(intervals, key=keyfn):
                yield i
    
class Tabix(IntervalFile, Index):
    def __init__(self, path, mode="bed"):
        if not os.path.exists(path + ".tbi"):
            if not path.endswith(".gz"):
                warnings.warn("""
                The file you are attempting to use as a 
                TabixGenomicIntervalIndex (%s) 
                is not BGZF compressed. It will be BGZF compressed and the 
                original file will be DELETED.""" % path)
            path = pysam.tabix_index(path, preset=mode)
        super(Tabix, self).__init__(path)
        self._tbi = pysam.Tabixfile(path)
   
    def fetch(self, query):
        return (Interval.parse_line(line, type="bed") for line in 
                self._tbi.fetch(query.contig, query.start, query.end))
                                                  
def open(path, indexed=False):
    """
    Generically open a genomic interval file as either an index or interval set.
    """
    if indexed:
        return Tabix(path)
    else:
        return IntervalFile(path)

@lru_cache()
def genome(name):
    import mysql.connector
    db = mysql.connector.connect(user="genome", password="",
                                 host="genome-mysql.cse.ucsc.edu",
                                 database=name)
    c = db.cursor()
    c.execute("SELECT chrom, size FROM chromInfo")
    return dict(c)

import argparse
import math
import operator
import os
import sys

import scipy.stats
import numpy
import pandas
import termcolor

def file_attribute(fn):
    """
    Decorator for a function that takes a file path as the first argument.
   
    Will memoize the results of the function and only recompute them if 
    the size or modification time of the file has changed.
    """
    import shelve
    raise NotImplementedError
    print(fn.__name__)
    print(fn.__module__)

#@file_attribute
#def foo(path):
#    print(path)
#    return hash(path)

# TODO:
# - What to do about contigs that are in query but not reference or
#   vice-versa? (currently ignored)
# - What to do about query points that occur before or after all
#   reference points? (currently ignored)

def overlap_enrichment_test(reference, query):
    raise NotImplementedError

def kolmogorov_smirnov_relative_distance_test(reference, query):
    reference = reference.midpoints()
    query = query.midpoints()
    distances = []
    for contig in set(reference.keys()) & set(query.keys()):
        q = query[contig]
        r = reference[contig]
        for i,ix in enumerate(numpy.searchsorted(r, q)):
            if ix in (0, len(r)):
                # FIXME?
                continue
            before,after = r[ix - 1], r[ix]
            d = min(abs(q[i] - before), abs(q[i] - after)) / abs(before - after)
            distances.append(d)
    distances = numpy.array(distances) * 2
    D, p = scipy.stats.kstest(distances, "uniform")
    sign = -1 if distances.mean() > 0.5 else 1
    return sign * -math.log10(p)

def ecdf_relative_distance_test(reference, query):
    raise NotImplementedError

def absolute_distance_test(reference, query):
    if True:
        raise NotImplementedError

    reference = reference.midpoints()
    query = query.midpoints()
    observed = []
    null = []

    for contig in set(reference.keys()) & set(query.keys()):
        q, r = query[contig], reference[contig]
        observed.extend(numpy.abs(r-e).min() for e in q)

    observed = numpy.array(observed)
    print(observed.mean())
   
def run(query_path, reference_paths, background=None):
   # FIXME: add multiprocessing
   columns = ["File", "KS Relative Distance [-log10(p)]"]
   is_p_value = [True]
   rows = []
   for reference_path in reference_paths:
      reference = dcaf.genome.open(reference_path)
      query = dcaf.genome.open(query_path)
      rows.append((reference_path,
                   kolmogorov_smirnov_relative_distance_test(reference, query)))
   df = pandas.DataFrame.from_records(rows, columns=columns, index="File")
   df.columns.is_p_value = is_p_value
   # FIXME: multiple hypothesis correction
   return df

def color_p_value(p_value, cutoff=0.05, log10=True, signed=True):
    """
    Given a p-value, return a ANSI terminal colored string if the 
    p-value is significant, or an ordinary string if not.
   
    By default, the input is expected to be a signed -log10 p-value, meaning
    that the number will be negative if there is underrepresentation and positive
    if overrepresentation.
   
    Colors:

    - green if significantly overrepresented
    - red if significantly underrepresented
    - yellow if unsigned and significant
    """
    sign = int(math.copysign(1.0, p_value))
    assert(sign in (-1,1))
    significant = abs(p_value) > pow(10,cutoff) if log10 else p_value < cutoff 
    significant = True
    if significant:
        if signed:
            color = "green" if sign == 1 else "red"
        else:
            color = "yellow"
        return termcolor.colored(p_value, color)
    else:
        return str(p_value)

@dcaf.util.entry_point
def genome_runner(args):
    """
    A command-line utility for statistically comparing sets of genomic intervals for 
    propensity to overlap for relative or absolute proximity to one another.

    For our previous publication on enrichment analysis:

    * http://www.ncbi.nlm.nih.gov/pubmed/22155868

    We also implemented several proximity analyses from the R GenometriCorr package:

    * http://www.ploscompbiol.org/article/info:doi/10.1371/journal.pcbi.1002529 (paper)
    * http://genometricorr.sourceforge.net/ (software package)
    """
    
    parser = argparse.ArgumentParser(description="GenomeRunner command-line version.")
    parser.add_argument("-b", "--background")
    parser.add_argument("-g", "--genome")
    parser.add_argument("query", nargs=1)
    parser.add_argument("references", nargs="+")
    opts = parser.parse_args(args)
   
    df = run(opts.query[0], opts.references)
    if os.fstat(0) == os.fstat(1):
        # Output is not redirected (is interactive); add coloring.
        for i, (column,is_p_value) in enumerate(zip(df.columns, df.columns.is_p_value)):
            if is_p_value:
                df[column] = list(map(color_p_value, df[column]))
    df.to_csv(sys.stdout, sep="\t")
