import codecs
import contextlib
import fcntl
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

def reader(handle):
    return codecs.getreader(locale.getpreferredencoding())(handle)

def writer(handle):
    return codecs.getwriter(locale.getpreferredencoding())(handle)

def _copy_file(src, dst):
    """
    Copy one text-mode file object to another.
    """
    try:
        for line in src:
            dst.write(line)
    except BrokenPipeError:
        pass

def download(url, text_mode=True, cache=True, cache_dir="/tmp/dcaf"):
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
    return gzip.open(path, mode)
 
def decode_stream(stream, encoding=locale.getpreferredencoding()):
    return codecs.getreader(encoding)(stream)

class Handle(dcaf.util.Proxy):
    """
    A wrapper for file-like objects (in text mode)
    that contains some additional convenience functions.
    """
    def __init__(self, handle, mode="rt"):
        if isinstance(handle, str):
            if handle.startswith("ftp://") or handle.startswith("http://"):
                handle = download(handle)
            elif handle.endswith(".gz"):
                handle = gzip.open(handle, mode)
            else:
                handle = open(handle, mode=mode)

        # FIXME: convert binary mode files to text-mode

        #if False:
        #    flags = fcntl.fcntl(handle.fileno(), fcntl.F_GETFL, 0)
        #    encoding = locale.getpreferredencoding()
        #    if os.O_RDONLY & flags:
        #        handle = codecs.getreader(encoding)(handle)
        #    else:
        #        handle = codecs.getwriter(encoding)(handle)
        
        super(Handle, self).__init__(handle)

    def __or__(self, cmd): 
        """
        Use the 'or' operator to pipe to subprocesses.
        """
        if isinstance(cmd, str):
            p = subprocess.Popen(cmd, bufsize=-1, shell=True,
                                 stdout=subprocess.PIPE, 
                                 stdin=subprocess.PIPE)
            o = writer(p.stdin)
            threading.Thread(target=_copy_file,
                             args=(self, o)).start()
            return Handle(reader(p.stdout))
        else:
            _copy_file(self, cmd)
    
    def __enter__(self, *args, **kwargs):
        return self._wrapped.__enter__(*args, **kwargs)

    def __exit__(self, *args, **kwargs):
        return self._wrapped.__exit__(*args, **kwargs)

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
