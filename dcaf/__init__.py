import os
import subprocess

__author__ = "Cory Giles and Mikhail Dozmorov"
__author_email__ = "mail@corygil.es"

# FIXME: This obviously won't work outside of a git clone
try:
    with open(os.devnull, "w") as devnull:
        __version__ = subprocess.check_output("git --git-dir %s describe --tags" % \
                                              os.path.join(os.path.dirname(__file__), "../.git"),
                                              stderr=devnull,
                                              shell=True).decode("ascii").strip()
except:
    try:
        from dcaf.version import __version__
    except ImportError:
        __version__ = "UNKNOWN"

_entry_points = {}

# Expose API modules

import dcaf.db
import dcaf.io
import dcaf.expression
