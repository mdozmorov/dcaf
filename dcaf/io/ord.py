from subprocess import Popen, PIPE
from collections import namedtuple

from . import ClosingMixin

Term = namedtuple("Term",
                  "accession name synonyms")

def _clean(text):
    return text.replace("\r","").replace("\n","").strip().strip('"')

class ORD(ClosingMixin):
    """
    Read terms and synonyms from the ORD MS Access database.
    
    Currently requires MDBTools to be installed (working on pure Python
    MDB driver).
    """
    ROW_DELIM = b"%%%%%"
    COL_DELIM = "$$$$$"

    def __init__(self, path):
        args = ["mdb-export", "-H",
                "-d", self.COL_DELIM, 
                "-R", self.ROW_DELIM,
                path, "tblObjectSynonyms"]
        self._pipe = Popen(args, stdout=PIPE)
        self._terms = {}
        
    def __iter__(self):
        p = self._pipe.stdout
        data = b""
        while True:
            new_data = p.read(4096)
            data += new_data
            splits = data.split(self.ROW_DELIM)
            for record in splits[:-1]:
                record = record.decode("utf8")
                fields = record.split(self.COL_DELIM)
                id = int(fields[0])
                accession = "ORD:" + str(id)
                name = _clean(fields[1])
                synonym = _clean(fields[3])
                if not id in self._terms:
                    self._terms[id] = Term(accession, name, [synonym])
                else:
                    self._terms[id].synonyms.append(synonym)
            if not new_data:
                break
            data = splits[-1]

        yield from self._terms.values()
                
    def close(self):
        self._pipe.stdout.close()

if __name__ == "__main__":
    import sys

    with ORD(sys.argv[1]) as h:
        for record in h:
            print(record)
