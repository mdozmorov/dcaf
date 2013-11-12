def test_go():    
    import dcaf.io

    url = "http://www.geneontology.org/ontology/obo_format_1_2/gene_ontology_ext.obo"
    with dcaf.io.download(url) as h:
        obo = dcaf.io.OBOFile(h)
        for term in obo:
            if term.id == "GO:0000001":
                assert term.name=="mitochondrion inheritance"
                assert "GO:0048308" in term.is_a
                assert "GO:0048311" in term.is_a
                return
