import pytest

import dcaf.io

@pytest.mark.parametrize("path,line_count", [
    (dcaf.io.data("manual_annotations_ursa.csv"), 14709),
    ("http://www.gnu.org/licenses/agpl-3.0.txt", 661)
])
def test_generic_open(path, line_count):
    import dcaf.io

    n = 0
    with dcaf.io.generic_open(path) as h:
        for line in h:
            n += 1
    assert n == line_count
