def test_medline_parser():
    import dcaf.io

    url = "http://www.nlm.nih.gov/databases/dtd/medsamp2013.xml"
    with dcaf.io.MedlineXMLFile(url) as h:
        n = 0
        for article in h:
            n += 1
    assert n == 156
