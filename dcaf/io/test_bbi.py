import pytest
import numpy.testing

import dcaf.io
import dcaf.io.bbi

@pytest.mark.parametrize("relative_path", [
    "test/bbi/test_bedGraph.bw",
    "test/bbi/test_varStep.bw"
])
def test_read_BigWig(relative_path):
    path = dcaf.io.data(relative_path)
    bw = dcaf.io.bbi.BigWigFile(path)

    # Basics
    assert len(bw) == 5430

    # Region query
    assert sum(1 for _ in bw.search("chr1", 10000, 11000)) == 21

    # Region summary
    summary = bw.summarize_region("chr1", 10000, 11000)
    
    assert summary.size == 1000
    assert summary.covered == 83
    numpy.testing.assert_allclose((summary.sum, summary.mean0, summary.mean),
                                  (-12.9303, -0.0129303, -0.155787),
                                  rtol=1e-5)

def test_read_BigBED():
    path = dcaf.io.data("test/bbi/ATF3.bb")
    bb = dcaf.io.bbi.BigBEDFile(path)

    # Basics
    assert len(bb) == 5459

    # Region search
    assert sum(1 for _ in bb.search("chr1", 0, 1000000)) == 3
    assert sum(1 for _ in bb.search("chrX", 0, 100000000)) == 39

def test_IntervalTree():
    raise NotImplementedError
#     itree = IntervalTree()
#     for i, (chrom, start, end) in enumerate(read_bed("scratch/knownGene.bed")):
#         itree.add(chrom,start,end,(chrom,start,end))
#     itree.build()
#     for item in itree.search("chr1", 10000,20000):
#         print(item)
#     for item in itree.search("chrX", 10000,200000):
#         print(item)
