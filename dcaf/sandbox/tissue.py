URSA_MEAN = 7.32
URSA_STDEV = 2.63

def infer_tissue(X):
    """
    Infer tissue for sample(s). The columns should be gene symbols.
    """
    if X.columns.dtype != 'O':
        raise Exception("The expression matrix columns are not object types, implying they do not contain gene symbols. Gene symbols are required for tissue imputation.")

    URSA_BIN = which("ursa")
    Xt = impute(X).T
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "24"
    with tempfile.NamedTemporaryFile() as h:
        X.T.dropna().to_csv(h.name, sep="\t")
        p = Popen([URSA_BIN, "-i", h.name, "-n", str(X.shape[0]), "-r", "1"],
                  cwd=dirname(URSA_BIN), stdout=PIPE, env=env)
        T = pandas.io.parsers.read_csv(p.stdout, skiprows=(6*X.shape[0] + 1),
                                       sep="\t", index_col=0, header=False,
                                       names=X.index)
    c.execute("""
    SELECT translate(term.name, ' ', '_'), term.id
    FROM term
    INNER JOIN ontology
    ON ontology.id=term.ontology_id
    WHERE ontology.namespace='BTO'""")
    m = dict([(k.replace("/", "--"),v) for k,v in c])
    T.index = [m[k] for k in T.index]
    T.index.name = "Term ID"
    T.columns.name = "Sample ID"
    return T
