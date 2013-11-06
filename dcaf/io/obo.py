"""
Simple reader for OBO (Open Biomedical Ontology) files.
"""

import collections

OBOTerm = collections.namedtuple("OBOTerm", [
    "id", "name", "synonym", "is_a", "part_of", "namespace"
])

class OBOFile(object):
    """
    A simple iterative reader for OBO (Open Biomedical Ontology) files.
    
    Once the OBOFile object is created by passing the __init__ function
    a file handle to the OBO file, the user can iterate over ontology terms
    using a for loop, etc.
    
    Currently, only a subset of the OBO specification is supported. 
    Namely, only the following attribute of [Term] entries:
    - id
    - name
    - synonym
    - is_a
    - part_of
    - namespace
    """

    def __init__(self, handle):
        """
        Create an OBOFile object.
        
        :param handle: File handle to the OBO file
        :type handle: A file handle in 'rt' mode
        """
        self._handle = handle
    
    def _make_term(self, attrs):
        if ("id" in attrs) and ("name" in attrs):
            for key in OBOTerm._fields:
                attrs.setdefault(key, [])
            return OBOTerm(**attrs)

    def __iter__(self):
        """
        Iterate over [Term] entries in the OBO file, consuming the 
        handle in the process.
        """
        is_term = False
        attrs = collections.defaultdict(list)
        for line in self._handle:
            line = line.strip()
            if line in ("[Term]", "[Typedef]", "[Instance]"):
                is_term = line == "[Term]"
                t = self._make_term(attrs)
                if t and is_term: yield t
                attrs = collections.defaultdict(list)
            elif line:
                try:
                    key, value = line.split(": ", 1)
                except ValueError:
                    continue

                if key == "id":
                    attrs["id"] = value
                elif key == "name":
                    attrs["name"] = value
                elif key == "is_a":
                    attrs["is_a"].append(value.split("!")[0].strip())
                elif key == "part_of":
                    attrs["part_of"].append(value.split("!")[0].strip())
                elif key == "namespace":
                    attrs["namespace"] = value
                elif key == "synonym":
                    value = value[1:]
                    attrs["synonym"].append(value[:value.find("\"")])
        t = self._make_term(attrs)
        if t and is_term: yield t
