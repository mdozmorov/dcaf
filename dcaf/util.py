"""
Utility functions.

.. moduleauthor:: Cory Giles <mail@corygil.es>
"""
import functools
import itertools
import os

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

class Proxy(object):
    """
    A class that wraps an object and forwards all attribute
    accesses to that object.

    Pointless on its own, this class is intended to be subclassed, and
    subclass methods will be called before the wrapped object's
    attributes. 
    """ 
    def __init__(self, wrapped): 
        self._wrapped = wrapped
        
    def __getattr__(self, attr):
        if attr in self.__dict__:
            return self.__dict___[attr]
        else:
            return getattr(self._wrapped, attr)
    
    def __hasattr__(self, attr):
        return attr in self.__dict__ or hasattr(self._wrapped)
    
    def __dir__(self):
        return list(set(self.__dict__.keys()) | set(dir(self._wrapped)))
  
def open_data(path):
    """
    Return a file handle to a data file (for data files included
    in the dcaf installation).

    :param path: Path to the data file, relative to the data root folder.
    :type path: str
    """
    return open(os.path.join(os.path.dirname(__file__), 
                             "..", "data", 
                             *os.path.split(path)))

def partition(n, seq):
    """
    Partition an iterator or collection into chunks of size n.
    
    .. warning::
        This function yields an *iterator* for each chunk,
        and as such, will not give correct results if each chunk
        is not completely realized before calling for the next chunk.
    """
    it = iter(seq)
    get_chunk = lambda: itertools.islice(it, n)
    chunk = get_chunk()
    while chunk:
        yield chunk
        chunk = get_chunk()
