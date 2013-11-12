def initialize_db():
    import dcaf.db

    db = DCAFConnection.from_configuration()
    db.load_schema()

    url = "http://www.geneontology.org/ontology/obo_format_1_2/gene_ontology_ext.obo"
    with dcaf.io.download(url) as h:
        db.import_obo(h, "GO", "Gene Ontology")
    
    with open("data/brenda.obo") as h:
        db.import_obo(h, "BTO", "Brenda Tissue Ontology")
