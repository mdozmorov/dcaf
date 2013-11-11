import codecs
import gzip
import hashlib
import locale
import os
import shutil
import urllib.request

import pandas

import dcaf

__all__ = ["download", "data", "read_matrix"]

def download(url, return_path=False,
             text_mode=True, cache=True, cache_dir="/tmp/dcaf"):
    """
    Download a URL, possibly caching and compressing the result.
    
    :param url: The URL to download
    :type url: str
    """
    # TODO: add cache expiration
    # TODO: make it tee during the initial download

    if not os.path.exists(cache_dir):
        os.makedirs(cache_dir)

    path = os.path.join(cache_dir, hashlib.md5(url.encode("ascii")).hexdigest())

    if not os.path.exists(path):
        i_handle = urllib.request.urlopen(url)
        if not cache:
            return io.TextIOWrapper(decode_stream(i_handle)) \
                if text_mode else i_handle
        o_handle = open(path, "wb") if url.endswith(".gz") else gzip.open(path, "wb")
        with i_handle:
            with o_handle:
                shutil.copyfileobj(i_handle, o_handle)

    mode = "rt" if text_mode else "rb"
    return path if return_path else gzip.open(path, mode)
 
def data(relative_path):
    """
    Return an absolute path to a data file (for data files included
    in the dcaf installation).

    :param path: Path to the data file, relative to the data root folder.
    :type path: str
    """        
    path = os.path.join(os.path.dirname(dcaf.__file__), 
                        "..", "data", 
                        *os.path.split(relative_path))
    return path if os.path.exists(path) else None

def read_matrix(path):
    """
    Read a tab-delimited matrix containing aligned row and column names 
    from a file path.
    
    :param path: Path to the matrix
    :type path: str
    :rtype: pandas.DataFrame
    """
    M = pandas.read_table(path, sep="\t", index_col=0)
    M.index = map(str, M.index)
    M.columns = map(str, M.columns)
    return M
