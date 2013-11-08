"""
A read-only API for BBI (BigWig and BigBED) files.

.. moduleauthor:: Cory Giles <mail@corygil.es>

.. todo::
    * add a "normalize" kwarg to "mean"
    * handle invalid BBI files
    * handle exons (in Python or C++?)
    * API for getting BBI header information
    * Parallelize the "same region across many files" case in BigWigSet
    * Finish BigBED class
    * Make BigWigSet read metadata
    * Make both BigBED and BigWig iterable
"""

import os

from collections import namedtuple

import numpy

cdef extern from "dcaf/mmap.hpp" namespace "dcaf":
    cdef cppclass MMFile:
        pass

cdef extern from "dcaf/bbi/bbi.hpp" namespace "dcaf::bbi":
    cdef struct BED:
        char* chrom
        int start, end
        char* name
        float score
        char strand
        char* rest

    cdef struct TotalSummary:
        long basesCovered
        double minVal, maxVal, sumData, sumSquares
    
    cdef struct Summary:
        int length, covered
        float sum, mean0, mean

    cdef cppclass BBIFile:
        TotalSummary* totalSummary

    cdef cppclass BigWigFile(BBIFile):
        BigWigFile(char* path)
        Summary summary(BED)
    
BBISummary = namedtuple("BBISummary",
                        ["basesCovered","minVal",
                         "maxVal","sumData","sumSquares"])

cdef class BBI:
    @property
    def basesCovered(self):
        return self.h.totalSummary.basesCovered

    @property
    def minVal(self):
        return self.h.totalSummary.minVal

    @property
    def maxVal(self):
        return self.h.totalSummary.maxVal

    @property
    def sumData(self):
        return self.h.totalSummary.sumData

    @property
    def sumSquares(self):
        return self.h.totalSummary.sumSquares

cdef class BigWig(BBI):
    """
    Query interface for BigWig files.
    """
    cdef BigWigFile *h

    def __cinit__(self, str path):
        bpath = path.encode("UTF-8")
        cdef const char* pPath = bpath
        self.h = new BigWigFile(pPath)

    def __dealloc__(self):
        del self.h

    def mean(self, str chrom, int start, int end):
        cdef BED bed
        chr = chrom.encode("UTF-8")
        bed.chrom = chr
        bed.start = start
        bed.end = end
        cdef Summary s = self.h.summary(bed)
        return s.mean0

class BigWigSet(object):
    """
    A higher-level interface for querying multiple BigWig files simultaneously.
    """
    def __init__(self, dir):
        self._handles = []
        self._mean_coverage = []
        dir = os.path.abspath(dir)
        for p in os.listdir(dir):
            if p.endswith(".bw"):
                h = BigWig(os.path.join(dir,p))
                self._handles.append(h)
                s = h.summary
                self._mean_coverage.append(s.sumData / s.basesCovered)
        self._mean_coverage = numpy.array(self._mean_coverage)
    
    def mean(self, str chrom, int start, int end, normalize=True):
        v = numpy.array([h.mean(chrom,start,end) for h in self._handles])
        if normalize:
            v /= self._mean_coverage
        return v
       
def test():
    bws = BigWigSet("/data/OMRF-Sjogren")
    print(bws._handles[0].sumData)
    #print bws.mean("chrX", 0, 100000)
