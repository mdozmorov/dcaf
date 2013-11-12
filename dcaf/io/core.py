import bz2
import codecs
import gzip
import hashlib
import locale
import os
import shutil
import urllib.request
import urllib.parse

import pandas

import dcaf

__all__ = ["generic_open", "download", 
           "data", "read_matrix", "ClosingMixin"]

class ClosingMixin(object):
    """
    Inherit from this mixin class to make a file-like
    object close upon exiting a 'with statement' block.

    Basically does the same thing as :py:func:`contextlib.closing`,
    except via inheritance.
    """

    def __enter__(self):
        return self

    def __exit__(self, *args):
        try:
            self.close()
        except:
            pass
    
    def close(self):
        raise NotImplementedError

def _is_url(path):
    """
    Guess whether this path is a URL.
    """
    parse = urllib.parse.urlparse(path)
    return all([parse.netloc,
                parse.scheme in ["http", "https", "ftp"]])

def generic_open(path, mode="rt"):
    """
    Open a file path, bzip2- or gzip-compressed file path, 
    or URL in the specified mode.

    Not all path types support all modes. For example, a URL is not
    considered to be writable. 
    
    :param path: Path or URL
    :type path: str
    :throws IOError: If the file cannot be opened in the given mode
    :throws FileNotFoundError: If the file cannot be found
    :rtype: :py:class:`io.IOBase` or :py:class:`io.TextIOBase`, 
      depending on the mode
    """
    # FIXME: detect zipped file based on magic number, not extension
    # FIXME: enable opening zipped file in text mode
    # FIXME: detect and unzip a zipped URL 
    # TODO: merge this and download (provide a cache_dir parameter or
    #   config setting)

    if _is_url(path):
        if "w" in mode:
            raise IOError("Cannot write to URL: '%s'")
        h = urllib.request.urlopen(path)
    else:
        h = open(path, mode)
    if path.endswith(".gz"):
        if not "b" in mode:
            raise IOError("Can't open zipped file in text mode (known bug).")
        h = gzip.GzipFile(fileobj=h)
    elif path.endswith(".bz2"):
        if not "b" in mode:
            raise IOError("Can't open zipped file in text mode (known bug).")
        h = bz2.BZ2File(h, mode=mode)
    return h

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
