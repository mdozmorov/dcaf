from dcaf.util import *

def test_which():
    assert(which("bash") in ("/usr/bin/bash", "/bin/bash"))
