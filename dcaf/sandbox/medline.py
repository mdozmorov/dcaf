def import_medline():
    import gzip
    import dcaf.io
    import lxml.etree

    xslt = lxml.etree.XSLT(lxml.etree.parse(dcaf.io.data("xslt/article.xsl")))
    root = xslt(dom)

    with gzip.open("/home/gilesc/data/text/medline13n0001.xml.gz") as h:
        dom = lxml.etree.parse(h)

    for row in io.StringIO(str(root)):
        print(row)
