def test_go():    
    import dcaf.io

    url = "http://www.geneontology.org/ontology/obo_format_1_2/gene_ontology_ext.obo"
    with dcaf.io.download(url) as h:
        obo = dcaf.io.OBOFile(h)
        for term in obo:
            if term.id == "GO:2001070":
                assert(term.name=="starch binding")
                assert("GO:0030247" in term.is_a)
                return
