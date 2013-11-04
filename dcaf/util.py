"""
Utility functions.

.. moduleauthor:: Cory Giles <mail@corygil.es>
"""

import functools
import sys

import dcaf

def entry_point(fn):
    """
    Register a function as a global entry point.
    
    Setuptools will create and install a wrapper script for this function. The
    function is assumed to take one argument, a list of command-line arguments.
    """
    @functools.wraps(fn)
    def wrapped(args=None):
        if args is None:
            args = sys.argv[1:]
        return fn(args)
    wrapped.__doc__ = fn.__name__ + "(argv)\n" + fn.__doc__
    dcaf._entry_points.append(wrapped)
    return wrapped
