"""
BigWig and BigBED readers in pure Python.

The detailed (byte-level) BigWig and BigBED specifications 
are described in the paper:

Kent WJ, Zweig AS, Barber G, Hinrichs AS, and Karolchik D.
BigWig and BigBed: enabling browsing of large distributed datasets.
Bioinformatics (2010) 26 (17): 2204-2207.

http://bioinformatics.oxfordjournals.org/content/26/17/2204.long
"""
import itertools
import mmap
import operator
import sys
import zlib

from collections import namedtuple
from functools import lru_cache
from struct import Struct, unpack

import numpy

from dcaf.genome import GenomicRegionRAMIndex

__all__ = ["BigWigFile", "BigBEDFile"]

# FIXME: don't assume native endianness
# TODO: (low priority) read from network

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

VarStepDatum = numpy.dtype([
    ("start", "i4"), ("value", "f4")
])

BEDGraph = namedtuple("BEDGraph", "chrom start end value")

# BigBED structure
BED = namedtuple("BED", "chrom start end name score strand rest")

# Region summary
RegionSummary = namedtuple("RegionSummary", "size covered sum mean0 mean") 

def read_struct(map, cls, count=None, offset=None):
    if offset is not None:
        map.seek(offset)

    _count = 1 if count is None else count
    size = cls._struct.size
    data = map.read(size * _count)

    items = [cls._make(cls._struct.unpack_from(data, offset=o)) \
             for o in range(0, size*_count, size)]
    return items[0] if count is None else items

class Tree(object):
    """
    Common superclass of BPTree and RTree, which share a similar structure.
    """
    def __init__(self, map, offset):
        self._map = map
        self._offset = offset + self.HeaderType._struct.size
        self.header = read_struct(map, self.HeaderType, offset=offset)
   
    def _read(self, offset):
        node = read_struct(self._map, Node, offset=offset)
        cls = self.LeafItem if node.isLeaf else self.NonLeafItem
        items = read_struct(self._map, cls, count=node.count)
        if node.isLeaf:
            yield from items
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

        # Read contigs from B+ Tree
        self._contig_by_id = {}
        self._contig_by_name = {}
        for leaf in BPTree(self._map, self.header.chromosomeTreeOffset):
            contig = Contig(leaf.key.decode("ascii").rstrip('\x00'), 
                            leaf.chromId, leaf.chromSize)
            self._contig_by_id[contig.id] = contig
            self._contig_by_name[contig.name] = contig

        # Read dataCount (# sections for BigWig, # elements in BigBED)
        self._map.seek(self.header.fullDataOffset)
        self.dataCount = int.from_bytes(self._map.read(4), byteorder="little")
        
    @property
    def _leaves(self):
        if not hasattr(self, "__leaves"):
            # Read index sections from RTree and cache them in an Interval Tree
            self.__leaves = GenomicRegionRAMIndex()
            for leaf in RTree(self._map, self.header.fullIndexOffset):
                for contig_id in range(leaf.startChromIx, leaf.endChromIx+1):
                    start = leaf.startBase if (contig_id == leaf.startChromIx) else 0
                    end = leaf.endBase if (contig_id == leaf.endChromIx) \
                          else self._contig_by_id[contig_id].size
                    self.__leaves.add(contig_id, start, end, leaf)
            self.__leaves.build()
        return self.__leaves

    def __del__(self):
        try:
            self._map.close()
            self._handle.close()
        except AttributeError:
            pass
    
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
        for leaf in self._leaves:
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
    """
    A BigWig format file. In this format, genomic regions 
    can be associated with a floating-point value. BigWig is a compressed, 
    indexed version of the BEDGraph or Wiggle formats.
    
    BigWig differs from BigBED in that genomic regions must be
    non-overlapping, each genomic region can only be associated with
    one value, and there is no additional information beyond the
    floating-point value associated with each region (such as name,
    strand, exons, etc.). 
    """ 

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
        data = data[BigWigSectionHeader._struct.size:]

        if section_header.type == BigWigFile.bedGraph:
            return numpy.frombuffer(data, dtype=BedGraphDatum)

        elif section_header.type == BigWigFile.varStep:
            elements = numpy.frombuffer(data, dtype=VarStepDatum)
            end = numpy.zeros(len(elements), dtype="f4")
            end[-1] = section_header.chromEnd
            end[:-1] = elements["start"][1:]
            columns = (elements["start"], end, elements["value"])
            records = numpy.core.records.fromarrays(columns,
                                                    dtype=BedGraphDatum)
            return records
        elif section_header.type == BigWigFile.fixedStep:
            msg = "FixedStep BigWig not currently supported."
            raise NotImplementedError(msg)
        
        else:
            raise IOError("BigWig data section type '%s' not a recognized value (corrupted or invalid file?)." % section_header.type)


    def _search_leaf_internal(self, leaf, contig_id, start, end):
        regions = self._read_section_internal(leaf)
        ix = (regions["start"] < end) & (start < regions["end"])
        return regions[ix]
    
    def _export_element(self, contig_id, element):
        return BEDGraph(self._contig_by_id[contig_id].name, *element)

    def _read_section(self, leaf):
        contig_id = leaf.startChromIx
        return (self._export_element(contig_id, e) \
                for e in self._read_section_internal(leaf))
    
    def _search_leaf(self, leaf, contig_id, start, end):
        return (self._export_element(contig_id, e) \
                for e in self._search_leaf_internal(leaf, contig_id, start, end))
   
    ######################################
    # Public API (BigWig specific methods)
    ######################################

    def _collapse_runs(self, elements):
        # The Kent utilities collapse adjacent elements with the
        # same value, so we do too for comparability
        prev = next(elements)
        for e in elements:
            if (prev.chrom == e.chrom) and \
               (prev.end == e.start) and \
               (abs(prev.value - e.value) < 1e-3):
                prev = BEDGraph(prev.chrom, prev.start, e.end, prev.value)
            else:
                yield prev
                prev = e
        yield prev
    
    def __iter__(self):
        it = super(BigWigFile, self).__iter__()
        return self._collapse_runs(it)
        
    def search(self, chrom, start, end):
        it = super(BigWigFile, self).search(chrom, start, end)
        return self._collapse_runs(it)

    def __len__(self):
        return sum(1 for _ in self)

    def summarize_region(self, chrom, start, end):
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

            if len(regions):
                # Truncate edge regions so that only the overlapping
                # segments are considered
                regions["start"][0] = max(regions["start"][0], start)
                regions["end"][-1] = min(regions["end"][-1], end)

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
    """
    A BigBED format file. This is essentially a BED file indexed by genomic
    region and containing various summary statistics. 
    """ 

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
    
    def __len__(self):
        return self.dataCount
