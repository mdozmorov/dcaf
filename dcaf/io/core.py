import codecs
import contextlib
import functools
import gzip
import hashlib
import locale
import os
import shutil
import subprocess
import sys
import threading
import urllib.request

import dcaf.util

def _is_binary(handle):
    """
    Attempt to infer whether a file handle is in binary mode.
    """
    return "b" in handle.mode if isinstance(handle.mode,str) else True

def _copy_file(src, dst):
    """
    Copy one file object to another, ignoring broken pipes and
    silently converting bytes to strings as necessary using the
    platform's default locale.
    """
    if _is_binary(dst) and not _is_binary(src):
        src = src.encode()
    elif not _is_binary(dst) and _is_binary(src):
        src = src.decode()
    try:
        shutil.copyfileobj(src, dst)
    except BrokenPipeError:
        pass

def _download(url, cache_dir="/tmp/dcaf"):
    """
    Download a URL, caching and compressing the result.
    
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
        o_handle = open(path, "wb") if url.endswith(".gz") else gzip.open(path, "wb")
        with i_handle:
            with o_handle:
                shutil.copyfileobj(i_handle, o_handle)

    return gzip.open(path, "r")
 
class Handle(dcaf.util.Proxy):
    """
    A wrapper for file-like objects that contains some additional
    convenience functions.
    """
    def __init__(self, handle, mode="r"):
        if isinstance(handle, str):
            if handle.startswith("ftp://") or handle.startswith("http://"):
                handle = _download(handle)
            elif handle.endswith(".gz"):
                handle = gzip.open(handle, mode)
            else:
                handle = open(handle, mode=mode)
        super(Handle, self).__init__(handle)

    def __or__(self, cmd): 
        """
        Use the 'or' operator to pipe to subprocesses.
        """
        if isinstance(cmd, str):
            p = subprocess.Popen(cmd, bufsize=-1, shell=True,
                                 stdout=subprocess.PIPE, 
                                 stdin=subprocess.PIPE)
            threading.Thread(target=_copy_file,
                             args=(self, p.stdin)).start()
            #return codecs.getreader("utf-8")(p.stdout)
            return Handle(p.stdout)
        else:
            _copy_file(self, cmd)
    
    def __iter__(self):
        return iter(self._wrapped)
    
    def __gt__(self, path):
        """
        Redirect this handle to a path.

        :param path: File path
        :type path: str
        """
        _copy_file(self, Handle(path, "wb"))
    
    def fields(self, sep="\t"):
        self.as_str()
        for line in self:
            yield line.strip().split(sep)

    def as_str(self):
        """
        Ensure this is a string stream. If it is a byte stream,
        decode it, otherwise do nothing.
        """
        if _is_binary(self):
            self.decode()
        return self
    
    def as_bytes(self):
        """
        Ensure this is a bytes stream. If it is a string stream,
        encode it, otherwise do nothing.
        """
        if not _is_binary(self):
            self.encode()
        return self

    def encode(self, encoding=locale.getpreferredencoding()):
        """
        Encode this stream using the given encoding.
        """
        self._wrapped = codecs.getwriter(encoding)(self._wrapped)
        return self
    
    def decode(self, encoding=locale.getpreferredencoding()):
        """
        Decode this stream using the given encoding.
        """
        self._wrapped = codecs.getreader(encoding)(self._wrapped)
        return self
