#!/usr/bin/env python

#
# Author:       Alex Reynolds
# Project:      Convert GFF3 to six-column UCSC BED
# Dependencies: Biopython (http://biopython.org/wiki/Biopython) and 
#               BCBio (https://github.com/chapmanb/bcbb)
#
# Example usage:
#
#  $ gff2bed.py < foo.gff | sort-bed -
#

import sys
from BCBio import GFF

def main(*args):
    max_target_lines = 1000
    requiredVersion = (2,5)
    checkInstallation( requiredVersion )

    gff_handle = sys.stdin
    for rec in GFF.parse(gff_handle, target_lines = max_target_lines):
        elem_chr = rec.id
        for feat in rec.features:
            elem_start = str(feat.location._start)
            end = feat.location._end
            elem_stop = str(feat.location._end.position - 1)
            elem_id = str(feat.id)
            try:
                elem_score = str(feat.qualifiers['score'][0])
            except KeyError:
                elem_score = '0'
            elem_strand = '+' if feat.strand == 1 else '-'
            print '\t'.join(["chr"+elem_chr, elem_start, elem_stop, elem_id, elem_score, elem_strand])

    return 0

def checkInstallation(rv):
    currentVersion = sys.version_info
    if currentVersion[0] == rv[0] and currentVersion[1] >= rv[1]:
        pass
    else:
        sys.stderr.write( "[%s] - Error: Your Python interpreter must be %d.%d or greater (within major version %d)\n" % (sys.argv[0], rv[0], rv[1], rv[0]) )
        sys.exit(-1)

    return 0

if __name__ == '__main__':
    sys.exit(main(*sys.argv))
