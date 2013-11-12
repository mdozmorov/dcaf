"""
An API for SQL databases containing various genomic data.

This module contains utilities to build and query a database
containing:

* Expression profiles from GEO 
* Text articles from PubMed/MEDLINE
* Metadata from Entrez Gene and RefSeq

Implementation details
~~~~~~~~~~~~~~~~~~~~~~

The core classes in this module are various subclasses of the
Connection class, which serve as thin wrappers around (actually,
proxies of) simple DB-API connections to various SQL databases.

Additionally, however, each Connection class implements extra methods
to query data from the given database and return results in an appropriate
formats; for example, expression datasets from the DCAFConnection will
return a pandas DataFrame, and datasets from UCSC will be returned as
genomic region objects. 

.. todo::
    * Add expression profiles from SRA
    * Various genomic feature tables from UCSC
    * Genomic coordinates for genes
"""

import io
import tarfile
import sys
import pandas

import psycopg2
import mysql.connector
from plumbum.cmd import tar, grep, cut, zcat, awk

import dcaf
import dcaf.ontology
from dcaf.util import Proxy, partition

# FIXME: Make DCAFConnection.expression id_type parameter an enum (Python 3.4 only)
# FIXME: Make PostgreSQLConnection and MySQLConnection

class Table(object):
    def __init__(self, db, name):
        self._db = db
        self._name = name
    
    def __iter__(self):
        return self._db("SELECT * FROM %s" % self.name)
        
    @property
    def name(self):
        return self._name
        
    @property
    def columns(self):
        rs = self._db("""
        SELECT column_name
        FROM information_schema.columns
        WHERE table_name = %s
        AND table_schema = 'public'
        ORDER BY ordinal_position""", self.name)
        return list(r[0] for r in rs)
    
    def copy(self, handle, **kwargs):
        """
        Copy rows from a file handle into this table 
        using the PostgreSQL COPY command.
        """
        c = self._db.cursor()
        c.copy_from(handle, self.name, **kwargs)
        self._db.commit()
       
    def truncate(self):
        """
        Delete all the rows in this table and dependent tables.
        """
        self._db("TRUNCATE %s CASCADE" % self.name)
    
class Connection(Proxy):
    """
    Wrapper for a database connection.
    """
    def __init__(self, user=None, password=None, 
                 schema=None, host=None, port=None):
        if self.DRIVER == "postgres":
            self._db = psycopg2.connect(user=user, password=password, 
                                        database=schema, 
                                        host=host, port=port)
        elif self.DRIVER == "mysql":
            self._db = mysql.connector.connect(user=user,
                                               password=password,
                                               host=host,
                                               port=port,
                                               database=schema)
        else:
            raise Exception("This database driver not supported.")

        self.schema = schema

        super(Connection,self).__init__(self._db)
    
    @classmethod
    def from_configuration(cls):
        """
        Connect to the database by searching for an appropriate 'dcaf.cfg' file.
        Any parameters from this file, if found, will override the 
        default connection parameters.
        """
        # FIXME: overridable kwargs
        cfg = dcaf.util.find_configuration()
        section = cls.CFG_SECTION
        return cls(user=cfg.get(section, "user"),
                   password=cfg.get(section, "password"),
                   schema=cfg.get(section, "schema"), 
                   host=cfg.get(section, "host"), 
                   port=cfg.get(section, "port"))

    def execute(self, sql, *params):
        """
        Execute SQL on this connection.
        
        :param sql: Raw SQL to be executed
        :type sql: str
        :param params: Parameters to be replaced in the format string
        :returns: A cursor object containing any applicable result set.
        """
        c = self._db.cursor()
        c.execute(sql, params)
        #self._db.commit()
        return c

    def __call__(self, sql, *params):
        """
        Shorthand for Connection.execute.
        """
        return self.execute(sql, *params)

    def execute_file(self, file):
        """
        Execute SQL contained in a file handle or file path.

        :param file: File path or file-like object
        :type file: str or handle
        :returns: A cursor object containing any applicable result set.
        """
        if not hasattr(file, "read"):
            file = open(file)
        with file:
            return self.execute(file.read())

    def table(self, name):
        return Table(self, name)
 
class DCAFConnection(Connection):
    """
    Wrapper for the dcaf PostgreSQL database.
    """
    DRIVER = "postgres"
    CFG_SECTION = "database"

    """
    Create a new connection to the dcaf database.

    If the connection is successful, but the schema 
    has not yet been loaded, then the DB schema will be loaded,
    and some small, core datasets will be imported.

    Large datasets like MEDLINE, GEO, etc, must be imported manually.
    """
    def load_schema(self):
        """
        (Re)load the database schema.
        """
        if not "taxon" in self.tables:
            self.execute_file(open(dcaf.io.data("sql/schema.sql")))
            self.execute_file(open(dcaf.io.data("sql/functions.sql")))

        self("""
        INSERT INTO ontology 
        VALUES (0, 'CORE', 'Core relations inherent to OBO format.')""")
        self("""
        INSERT INTO term (id, ontology_id, name) 
        VALUES (0, 0, 'is_a')""")
        self("""INSERT INTO term (id, ontology_id, name) 
        VALUES (1, 0, 'part_of')""")

        self.commit()

    def initialize(self):
        """
        Load the schema, then download and import core datasets.
        """
        self.load_schema()

        NCBI_BASE = "ftp://ftp.ncbi.nlm.nih.gov"
        
        # Load NCBI taxa into 'taxon' table
        url = NCBI_BASE + "/pub/taxonomy/taxdump.tar.gz"
        path = dcaf.io.download(url, return_path=True)
        self.table("taxon").copy(
            (tar["Oxzf", path, "names.dmp"] | grep["scientific"] \
             | cut["-f1,3"]).popen().stdout)

        # Load Entrez Gene IDs, names, symbols, etc, into 'gene'
        path = dcaf.io.download(NCBI_BASE + "/gene/DATA/gene_info.gz", 
                                return_path=True)
        self.table("gene").copy(
            (zcat[path] | grep["-v", "^#"] | cut["-f1-3,12"] \
            | awk["$1 != 366646"]).popen().stdout)

    def import_medline(self, folder):
        """
        Import MEDLINE XML files into the database.
        """
        for file in os.listdir(folder):
            if file.endswith(".xml.gz"):
                path = os.path.join(folder, file)

    def import_obo(self, handle, namespace, description):
        """
        Import an Open Biomedical Ontology (OBO) file into the database.
        """
        # is_a = 0
        # part_of = 1
        ontology_id = next(self("""
          INSERT INTO ontology (namespace, description) 
          VALUES (%s,%s) 
          RETURNING id""", namespace, description))[0]

        by_accession = {}
        relations = []
        obo = dcaf.io.OBOFile(h)
        for term in obo:
            term_id = next(self("""
              INSERT INTO term (ontology_id, accession, name) 
              VALUES (%s,%s,%s) RETURNING id""",
                                ontology_id, term.id, term.name))[0]
            by_accession[term.id] = term_id
            
            for parent in term.is_a:
                relations.append((term.id, parent, 0))
            for parent in term.part_of:
                relations.append((term.id, parent, 1))

            for synonym in term.synonym:
                self("""
                INSERT INTO synonym (term_id, synonym)
                VALUES (%s,%s)""", term_id, synonym)

        self.commit()
        
        for (child, parent, type) in relations:
            child = by_accession[child]
            parent = by_accession[parent]
            if child and parent:
                self("""
                  INSERT INTO term_relation (agent, target, relation)
                  VALUES (%s,%s,%s)""", child, parent, type)

        self.commit()

    @property
    def tables(self):
        q = """
        SELECT table_name 
        FROM information_schema.tables
        WHERE table_schema = 'public'
        """
        return list(row[0] for row in self.execute(q))
        
  
    def ontology(self, taxon_id=9606, namespace="GO"):
        """
        Load an ontology from this database into a directed graph.
        """
        return dcaf.ontology.Ontology(self._db.cursor(), 
                                      taxon_id=taxon_id, namespace=namespace)
        
    def genes(self, taxon_id):
        c = self("""
        SELECT id, symbol, name 
        FROM gene 
        WHERE taxon_id=%s 
        ORDER BY id""", (taxon_id,))
        return pandas.DataFrame.from_items([(row[0], row) for row in c], 
                                           columns=["id", "symbol", "name"],
                                           orient="index")

    def expression(self, taxon_id=9606, id_type="entrez", require_age=True, 
                   require_gender=False, shuffle=False, limit=None):
        """
        Fetch an expression set corresponding to the given criteria.
        """
        assert(id_type in ("entrez", "symbol"))

        genes = self.genes(taxon_id)
        query = """
        SELECT sample.id, sample.age, sample.gender, expression.data 
        FROM expression 
        INNER JOIN sample 
        ON expression.sample_id=sample.id 
        INNER JOIN platform
        ON sample.platform_id=platform.id
        WHERE platform.taxon_id=%s"""
        if require_age:
            query += "\nAND sample.age IS NOT NULL"
        if require_gender:
            query += "\nAND sample.gender IS NOT NULL"
        if limit:
            query += "\tLIMIT " + str(limit)
        if shuffle:
            query = "SELECT * FROM (%s) AS subq ORDER BY random()" % query

        c = self._db.cursor("expression")
        c.itersize = 1
        c.execute(query, (taxon_id,))
        samples, age, gender, expression = zip(*c)
        X = pandas.DataFrame.from_records(list(expression), 
                                          index=samples, columns=genes.index)
        X.index.name = "Sample ID"
        X.columns.name = "Gene ID"
        P = pandas.DataFrame({"age": age, "gender": gender}, index=samples)
        P.index.name = "Sample"

        if id_type == "symbol":
            X = self._collapse_to_symbols(X, taxon_id)

        return X, P

    def sample_expression(self, accession):
        """
        Get expression for a single sample.
        """
        c = self("""
        SELECT data FROM expression
        INNER JOIN sample
        ON sample.id=expression.sample_id
        WHERE sample.accession=%s""", (accession,))
        try:
            return Series(c.__next__()[0], 
                          index=self.genes(9606).index)
        except StopIteration:
            return None

    def _symbol_map(self, taxon_id):
        m = {}
        return dict(self("""
        SELECT id, symbol 
        FROM gene 
        WHERE taxon_id=%s AND symbol IS NOT NULL""", taxon_id))

    def _collapse_to_symbols(self, X, taxon_id, axis=1, min_pct=0.2):
        # FIXME: collapse by max mean instead of max
        thresh = int(X.shape[0] * min_pct)
        return X.groupby(self._symbol_map(taxon_id), 
                         axis=1).max().dropna(thresh=thresh, axis=1)

class UCSCConnection(Connection):
    """
    A connection to the (MySQL) database underlying the UCSC Genome Browser.
    """
    DRIVER = "mysql"
    CFG_SECTION = "UCSC"
