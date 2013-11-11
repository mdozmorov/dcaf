"""
BigWig and BigBED readers in pure Python.
"""
import mmap
import sys
import zlib

from collections import namedtuple
from functools import lru_cache
from struct import Struct, unpack, unpack_from, calcsize

import numpy

# FIXME: don't assume native endianness

class IntervalNode(object):
    def __init__(self, intervals):
        midpoint = int(len(intervals) / 2)
        start, end, data = intervals[midpoint]
        self._start = start
        self._end = end
        self._data = data
        self._left = None
        self._right = None
        self._max_end = self._end
        if midpoint > 0:
            self._left = IntervalNode(intervals[:midpoint])
            self._max_end = max(self._end, self._left._max_end)
        if midpoint < len(intervals) - 1:
            self._right = IntervalNode(intervals[midpoint+1:])
            self._max_end = max(self._end, self._right._max_end)

    def search(self, start, end):
        if start > self._max_end:
            return
        if self._left:
            yield from self._left.search(start, end)
        if (start < self._end) and (end > self._start):
            yield self._data
        if end < self._start:
            return
        if self._right:
            yield from self._right.search(start, end)

class IntervalTree(object):
    def __init__(self):
        self._intervals = {}
        self._roots = {}
        self._built = False

    def add(self, chrom, start, end, data=None):
        if self._built:
            raise Exception("Can't currently mutate a constructed IntervalTree.")
        self._intervals.setdefault(chrom, [])
        self._intervals[chrom].append((start, end, data))
    
    def build(self):
        assert(self._intervals)
        for chrom, intervals in self._intervals.items():
            intervals.sort()
            self._roots[chrom] = IntervalNode(intervals)
        self._built = True
    
    def search(self, chrom, start, end):
        if not self._built:
            raise Exception("Must call IntervalTree.build() before using search()")
        return self._roots[chrom].search(start, end)
 
def defstruct(name, format, fields):
    cls = namedtuple(name, fields)
    cls._struct = Struct(format)
    return cls

# General BBI file structures

BBIHeader = defstruct("BBIHeader", "=IHHQQQHHQQIQ",
                      "magic version zoomLevels \
                      chromosomeTreeOffset \
                      fullDataOffset fullIndexOffset \
                      fieldCount definedFieldCount \
                      autoSqlOffset totalSummaryOffset \
                      uncompressBufSize reserved")

TotalSummary = defstruct("TotalSummary", "=Qdddd",
                         "basesCovered minVal \
                         maxVal sumData sumSquares")

Contig = namedtuple("Contig", "name id size")

# Common Node metadata for BPTree and RTree
Node = defstruct("Node", "=??H", "isLeaf reserved count")

# BPTree structures
BPTreeHeader = defstruct("BPTreeHeader", "=IIIIQQ",
                         "magic blockSize keySize \
                         valSize itemCount reserved")

# RTree structures
RTreeHeader = defstruct("RTreeHeader", "=IIQIIIIQII",
                        "magic blockSize itemCount \
                        startChromIx startBase \
                        endChromIx endBase endFileOffset \
                        itemsPerSlot reserved")

RNonLeafItem = defstruct("RNonLeafItem", "=IIIIQ",
                         "startChromIx startBase \
                         endChromIx endBase childOffset")

RLeafItem = defstruct("RLeafItem", "=IIIIQQ",
                      "startChromIx startBase \
                      endChromIx endBase \
                      dataOffset dataSize")

# BigWig structures
BigWigSectionHeader = defstruct("BigWigSectionHeader",
                                "=IIIIIBBH",
                                "chromId chromStart \
                                chromEnd itemStep itemSpan \
                                type reserved itemCount")

BedGraphDatum = numpy.dtype([
    ("start", "i4"), ("end", "i4"), ("value", "f4")
])

BEDGraph = namedtuple("BEDGraph", "chrom start end value")

# BigBED structure
BED = namedtuple("BED", "chrom start end name score strand rest")

# Region summary
RegionSummary = namedtuple("RegionSummary", "size covered sum mean0 mean") 

def read_struct(map, cls, count=1, offset=None):
    if offset is not None:
        map.seek(offset)

    size = cls._struct.size
    data = map.read(size * count)

    items = [cls._make(cls._struct.unpack_from(data, offset=o)) \
             for o in range(0, size*count, size)]
    return items[0] if count == 1 else items

class Tree(object):
    """
    Common superclass of BPTree and RTree, which share a similar structure.
    """
    def __init__(self, map, offset):
        self._map = map
        self._offset = offset + self.HeaderType._struct.size
        self.header = read_struct(map, self.HeaderType, offset=offset)
   
    def _read(self, offset):
        self._map.seek(offset)
        node = read_struct(self._map, Node)
        cls = self.LeafItem if node.isLeaf else self.NonLeafItem
        items = read_struct(self._map, cls, count=node.count)
        if node.isLeaf:
            yield from iter(items)
        else:
            for item in items:
                yield from self._read(item.childOffset)
    
    def __iter__(self):
        yield from self._read(self._offset)
    
class BPTree(Tree):
    """
    The chromosome index, which is stored in a B-Plus Tree.
    """
    HeaderType = BPTreeHeader

    def __init__(self, map, offset):
        super(BPTree, self).__init__(map, offset)
        assert(self.header.magic==0x78CA8C91)
        self.LeafItem = defstruct("LeafItem", "=%ssII" % self.header.keySize, 
                                  "key chromId chromSize")
        self.NonLeafItem = defstruct("NonLeafItem", "=%ssQ" % self.header.keySize, 
                                     "key childOffset")

class RTree(Tree):
    """
    The interval index, which is stored in an R Tree.
    """
    HeaderType = RTreeHeader
    LeafItem = RLeafItem
    NonLeafItem = RNonLeafItem
    
    def __init__(self, map, offset):
        super(RTree, self).__init__(map, offset)
    
class ContigNotFound(Exception):
    def __init__(self, name, path):
        msg = "Contig '%s' not found in '%s'" % (name, path)
        super(ContigNotFound, self).__init__(msg)

class BBIFile(object):
    """
    Common superclass of BigWig and BigBed files.
    """
    def __init__(self, path):
        self._path = path
        self._handle = open(path, "r+b")
        self._map = mmap.mmap(self._handle.fileno(), 0)
        self.header = read_struct(self._map, BBIHeader)
        self.summary = read_struct(self._map, TotalSummary,
                                   offset=self.header.totalSummaryOffset)

        self._contig_by_id = {}
        self._contig_by_name = {}
        for leaf in BPTree(self._map, self.header.chromosomeTreeOffset):
            contig = Contig(leaf.key.decode("ascii").rstrip('\x00'), 
                            leaf.chromId, leaf.chromSize)
            self._contig_by_id[contig.id] = contig
            self._contig_by_name[contig.name] = contig

        self._leaves = IntervalTree()
        for leaf in RTree(self._map, self.header.fullIndexOffset):
            for contig_id in range(leaf.startChromIx, leaf.endChromIx+1):
                start = leaf.startBase if (contig_id == leaf.startChromIx) else 0
                end = leaf.endBase if (contig_id == leaf.endChromIx) else self._contig_by_id[contig_id].size
                self._leaves.add(contig_id, start, end, leaf)
        self._leaves.build()
                                        
    def __del__(self):
        self._map.close()
        self._handle.close()
    
    ########################
    # Decoding data sections
    ########################

    def _decompress_section(self, leaf):
        self._map.seek(leaf.dataOffset)
        data = self._map.read(leaf.dataSize)
        if self.header.uncompressBufSize > 0:
            data = zlib.decompress(data, 15, self.header.uncompressBufSize)
        return data

    def _read_section(self, leaf):
        """
        Uncompress the data under this leaf, and read the data contained therein.
        """
        data = self._decompress_section(leaf)
        return self._read_section_data(leaf, data)
    
    def _read_section_data(self, leaf, data):
        raise NotImplementedError
    
    ###################
    # Utility functions
    ###################

    def _get_contig_id(self, name):
        contig = self._contig_by_name.get(name)
        if contig is None:
            raise ContigNotFound(name, self._path)
        return contig.id

    #######################
    # Searching for regions
    #######################

    def _search_index(self, contig_id, start, end):
        """
        Search the RTree for RLeafItems overlapping the given region.
        """
        yield from self._leaves.search(contig_id, start, end)
    
    def _search(self, contig_id, start, end):
        for leaf in self._search_index(contig_id, start, end):
            yield from self._search_leaf(leaf, contig_id, start, end)
    
    def _search_leaf(self, leaf, contig_id, start, end):
        # Note, every section in BigWig files contains all the same
        # contig ID, whereas BigBED can have mixed contig IDs.
        # Thus, the _search_leaf method is implemented by subclasses
        # so that BigWig files can be searched more efficiently.
        raise NotImplementedError

    ############
    # Public API
    ############
    def __iter__(self):
        """
        Iterate through the elements (BED or BEDGraph) in this BBI file.
        """
        for leaf in self._rtree:
            yield from self._read_section(leaf)
    
    def search(self, contig, start, end):
        """
        Return the elements (BED or BEDGraph) in this BBI file
        that overlap the query region.
        """
        id = self._get_contig_id(contig)
        yield from self._search(id, start, end)

    @property
    def contigs(self):
        """
        Return a list of the contigs in this BBI file.
        """
        return list(self._contig_by_id.values())

    def __repr__(self):
        return "<%s: %s>" % (str(type(self)), self._path)
 

class BigWigFile(BBIFile):
    # Data section formats
    bedGraph = 1
    varStep = 2
    fixedStep = 3

    def __init__(self, path):
        super(BigWigFile, self).__init__(path)
        assert(self.header.magic==0x888FFC26)
        cache = lru_cache(maxsize=16)
        self._read_section_internal = cache(self._read_section_internal)
    
    def _read_section_internal(self, leaf):
        data = self._decompress_section(leaf)
        section_header_values = BigWigSectionHeader._struct.unpack_from(data)
        section_header = BigWigSectionHeader._make(section_header_values)
        assert(section_header.type in (1,2,3))
        data = data[BigWigSectionHeader._struct.size:]

        if section_header.type == BigWigFile.bedGraph:
            return numpy.frombuffer(data, dtype=BedGraphDatum)
        else:
            msg = "Only BedGraph type BigWig currently supported."
            raise NotImplementedError(msg)

    def _search_leaf_internal(self, leaf, contig_id, start, end):
        regions = self._read_section_internal(leaf)
        ix = (regions["start"] >= start) & (regions["end"] <= end)
        return regions[ix]
    
    def _export_element(self, contig_id, element):
        return BEDGraph(self._contig_by_id[contig_id].name, *element)

    def _read_section(self, leaf):
        contig_id = leaf.startChromIx
        return (self._export_element(contig_id, e) \
                for e in self._read_section_internal(leaf))
    
    def _search_leaf(self, leaf, contig_id, start, end):
        return (self._export_element(contig_id, e) \
                for e in self._search_leaf_internal(leaf, data))
   
    ######################################
    # Public API (BigWig specific methods)
    ######################################

    def summarize_region(self, chrom, start, end):
        # FIXME: doesn't properly handle partially overlapping regions
        length = end - start
        assert(length > 0)

        try:
            contig_id = self._get_contig_id(chrom)
        except ContigNotFound:
            return RegionSummary(length, 0, 0, 0, 0)

        bases_covered = 0
        sum_values = 0
        for leaf in self._search_index(contig_id, start, end):
            regions = self._search_leaf_internal(leaf, contig_id, start, end)
            lengths = regions["end"] - regions["start"]
            sum_values += (lengths * regions["value"]).sum()
            bases_covered += lengths.sum()
        if sum_values == 0:
            mean0, mean = 0, 0
        else:
            mean0 = sum_values / length
            mean = sum_values / bases_covered
        return RegionSummary(length, bases_covered, sum_values, mean0, mean)

class BigBEDFile(BBIFile):
    def __init__(self, path):
        super(BigBEDFile, self).__init__(path)
        assert(self.header.magic==0x8789F2EB)
    
    def _search_leaf(self, leaf, contig_id, start, end):
        contig_name = self._contig_by_id[contig_id].name
        for bed in self._read_section(leaf):
            if contig_name == bed.chrom:
                if (bed.start < end) and (start < bed.end):
                    yield bed

    def _read_section_data(self, leaf, data):
        records = []
        while data:
            contig, start, end = unpack("=III", data[:12])
            data = data[12:]
            _null = data.find(b'\0')
            fields = data[:_null].split(b'\t')
            name = fields[0] if fields else sys.intern('.')
            score = float(fields[1]) if len(fields) > 1 else numpy.nan
            strand = fields[2] if len(fields) > 2 else sys.intern('.')
            rest = tuple(fields[3:]) if len(fields) > 3 else None
            record = BED(self._contig_by_id[contig].name, 
                         start, end, name, score, strand, rest)
            records.append(record)
            data = data[(_null+1):]
        return records

def read_bed(path):
    h = open(path)
    for line in h:
        chrom, start, end = line.split("\t")[:3]
        yield chrom, int(start), int(end)
    h.close()

def testBigWig():
    bw = BigWigFile("/home/gilesc/data/RNAseq/bw/DRR001622.bw")
    with open("scratch/knownGene.bed") as h:
        for i,line in enumerate(h):
            chrom, start, end = line.split("\t")[:3]
            print(chrom, start, end, 
                  bw.summarize_region(chrom, int(start), int(end)).mean)
            #if i == 10000:
            #    break
    #print(bw.summarize_region("chr1", 0, 100000))
    #print(bw.summary)
    #data = bw._read_section(next(iter(bw._rtree)))
    #print(data["value"])

def testBigBED():
    bb = BigBEDFile("/home/gilesc/labcode/cbio/test/ATF3.bb")
    for bed in bb.search("chr1", 0, 100000000):
        print(bed)

def testIntervalTree():
    itree = IntervalTree()
    for i, (chrom, start, end) in enumerate(read_bed("scratch/knownGene.bed")):
        itree.add(chrom,start,end,(chrom,start,end))
    itree.build()
    for item in itree.search("chr1", 10000,20000):
        print(item)
    for item in itree.search("chrX", 10000,200000):
        print(item)

testBigWig()
#testBigBED()
#testIntervalTree()

#import dcaf.io.bbi
#bws = dcaf.io.bbi.BigWigSet("~/data/RNAseq/bw/")
#print(bws._handles[0].sumData)
