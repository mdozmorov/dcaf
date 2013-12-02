"""
Parser for MSigDB XML format.
"""

import sys
from collections import namedtuple
from xml.etree import ElementTree

from dcaf.io import generic_open as open, ClosingMixin

__all__ = ["MSigDB", "GeneSet"]

GeneSet = namedtuple("GeneSet", "name description long_description accession organism category contributor genes gene_symbols")

class MSigDB(ClosingMixin):
    def __init__(self, path):
        self._handle = open(path)
        self._tree = ElementTree.iterparse(self._handle)
        
    def close(self):
        self._handle.close()
    
    def __iter__(self):
        for event, elem in self._tree:
            if not elem.tag == "GENESET":
                continue
            organism = elem.get("ORGANISM")
            genes = elem.get("MEMBERS_EZID")
            if not genes:
                continue
            genes = list(set(map(int, genes.split(","))))
            yield GeneSet(elem.get("STANDARD_NAME"),
                          elem.get("DESCRIPTION_BRIEF"),
                          elem.get("DESCRIPTION_FULL"),
                          elem.get("SYSTEMATIC_NAME"),
                          organism,
                          elem.get("CATEGORY_CODE"),
                          elem.get("CONTRIBUTOR"),
                          genes,
                          elem.get("MEMBERS_SYMBOLIZED", "").split(","))
        self.close()

if __name__ == "__main__":
    with MSigDB(sys.argv[1]) as h:
        for item in h:
            print(item.name, item.description)

