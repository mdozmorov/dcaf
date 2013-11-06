import codecs
import gzip
import hashlib
import locale
import os
import shutil
import urllib.request

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
            return decode_stream(i_handle) if text_mode else i_handle
        o_handle = open(path, "wb") if url.endswith(".gz") else gzip.open(path, "wb")
        with i_handle:
            with o_handle:
                shutil.copyfileobj(i_handle, o_handle)

    mode = "rt" if text_mode else "rb"
    return path if return_path else gzip.open(path, mode)
 
def decode_stream(stream, encoding=locale.getpreferredencoding()):
    return codecs.getreader(encoding)(stream)
