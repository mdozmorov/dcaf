"""
Utility functions.

.. moduleauthor:: Cory Giles <mail@corygil.es>
"""
import configparser
import functools
import itertools
import locale
import os
import subprocess
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

def which(exe):
    """
    Find the full path to an executable (calls the shell command).
    """
    encoding = locale.getpreferredencoding()
    path = subprocess.check_output(["which", exe]).decode(encoding).strip()
    return os.path.realpath(path)

def multimap(pairs):
    """
    Given an iterable of pairs, construct a dict mapping the first elements
    to sets of second elements.
    """
    m = {}
    for k,v in pairs:
        if not k in m:
            m[k] = set()
        m[k].add(v)
    return m

def coo_to_df(triplets):
    """
    Create a SparseDataFrame from a sequence of (row,col,value) triplets.
    """
    data = defaultdict(dict)
    for row, col, val in triplets:
        data[col][row] = val
    return SparseDataFrame(data)

class ConfigurationNotFound(Exception):
    def __init__(self):
        super(ConfigurationNotFound, self).__init__()

def find_configuration():
    """
    Search sys.path for a 'dcaf.cfg' file. If found, return an open ConfigParser
    for it. Otherwise, return a ConfigParser with the default configuration.
    """
    parser = configparser.ConfigParser(allow_no_value=True)
    with open_data("defaults.cfg") as h:
        parser.read_file(h)

    for folder in sys.path:
        if os.path.exists(folder) and os.path.isdir(folder):
            for file in os.listdir(folder):
                if file == "dcaf.cfg":
                    path = os.path.join(folder, file)
                    parser.read(path)
                    return parser
    return parser
