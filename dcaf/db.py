"""
An API for SQL databases containing various genomic data.

The core of this module are various subclasses of the Connection class,
which serve as thin wrappers around (actually, proxies of) simple DB-API
connections to various SQL databases. 

Additionally, however, each Connection class implements extra methods
to query data from the given database and return results in an appropriate
formats; for example, expression datasets from the DCAFConnection will
return a pandas DataFrame, and datasets from UCSC will be returned as
genomic region objects. 
"""
import io
import tarfile

import psycopg2
import mysql.connector

import dcaf
import dcaf.ontology
from dcaf.util import Proxy, partition
from dcaf.io import Handle

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
        """
        for i,group in enumerate(partition(100000, handle)):
            print(i)

            buffer = io.BytesIO()
            for line in group:
                buffer.write(line)
            buffer.seek(0)

            c.copy_from(buffer, self.name, **kwargs)
            self._db.commit()
        """
        
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
        self._db.commit()
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
            self.execute_file(dcaf.util.open_data("sql/schema.sql"))
            self.execute_file(dcaf.util.open_data("sql/functions.sql"))

    def _import_rows(self, table_name, rows):
        table = self.table(table_name)
        table.truncate()
        sql = "INSERT INTO %s VALUES (%s)" % \
              (table_name, ",".join(["%s" for _ in table.columns]))

        c = self.cursor()
        for group in dcaf.util.partition(1000000, rows):
            c.executemany(sql, group)
        self.commit()

    def initialize(self):
        """
        Load the schema, then download and import core datasets.
        """
        self.load_schema()
        BASE = "ftp://ftp.ncbi.nlm.nih.gov"

        print("Loading taxon ...")
        def taxon_rows():
            url = BASE + "/pub/taxonomy/taxdump.tar.gz"
            with dcaf.io.download(url, text_mode=False) as stream:
                tar = tarfile.TarFile(fileobj=stream)
                with tar.extractfile("names.dmp") as h:
                    for line in dcaf.io.decode_stream(h):
                        if "scientific" in line:
                            fields = line.split("\t")
                            yield fields[0], fields[2]

        self._import_rows("taxon", taxon_rows())

        print("Loading gene ...")
        def gene_rows():
            url = BASE + "/gene/DATA/gene_info.gz"
            with dcaf.io.download(url) as h:
                h.__next__()
                for line in h:
                    line = line.strip()
                    if not line:
                        continue
                    fields = line.split("\t")
                    row = (item if item != "-" else None
                           for item in fields[:3] + [fields[11]])
                    if row[0] != "366646":
                        yield row
        
        self._import_rows("gene", gene_rows())
       
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
        return DataFrame.from_items([(row[0], row) for row in c], 
                                    columns=["id", "symbol", "name"],
                                    orient="index")

    def expression(self, 
                   taxon_id=9606, 
                   id_type="entrez",
                   require_age=True, 
                   require_gender=False, 
                   limit=None):
        """
        Get an expression set corresponding to the given criteria.
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

        c = self._db.cursor("expression")
        c.itersize = 1
        c.execute(query, (taxon_id,))
        samples, age, gender, expression = zip(*c)
        X = DataFrame.from_records(list(expression), 
                                   index=samples, columns=genes.index)
        X.index.name = "Sample ID"
        X.columns.name = "Gene ID"
        P = DataFrame({"age": age, "gender": gender}, index=samples)
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
    DRIVER = "mysql"
    CFG_SECTION = "UCSC"

if __name__ == "__main__":
    db = DCAFConnection.from_configuration()
    db.initialize()
