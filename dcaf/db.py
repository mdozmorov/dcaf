"""
API for SQL database containing various genomic data.
"""
import io
import tarfile

import psycopg2

from dcaf.util import Proxy, partition
from dcaf.io import Handle

class Table(object):
    def __init__(self, db, name):
        self._db = db
        self._name = name
    
    def __iter__(self):
        return self._db("SELECT * FROM %s" % self.name)
        
    @property
    def name(self):
        return self._name
        
    def columns(self):
        raise NotImplementedError
    
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
    Wrapper for a connection to the dcaf PostgreSQL database.
    """
    def __init__(self, user=None, password=None, 
                 database=None, host=None, port=None):
        """
        Create a new connection to the dcaf database.
        
        If the connection is successful, but the schema 
        has not yet been loaded, then the DB schema will be loaded,
        and some small, core datasets will be imported.
        
        Large datasets like MEDLINE, GEO, etc, must be imported manually.
        """
        self._db = psycopg2.connect(user=user, password=password, 
                                    database=database, 
                                    host=host, port=port)
        super(Connection,self).__init__(self._db)

        if not "taxon" in self.tables:
            self.initialize()

    def initialize(self):
        """
        Load the database schema and core datasets.
        """
        self.execute_file(dcaf.util.open_data("sql/schema.sql"))
        self.execute_file(dcaf.util.open_data("sql/functions.sql"))
        
    @property
    def tables(self):
        q = """
        SELECT table_name 
        FROM information_schema.tables
        WHERE table_schema = 'public'
        """
        return list(row[0] for row in self.execute(q))
        
    def table(self, name):
        return Table(self, name)
        
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

if __name__ == "__main__":
    import sys
    BASE = "ftp://ftp.ncbi.nlm.nih.gov"

    db = Connection(database="testdb")

    # print("Loading taxon ...")
    # db.table("taxon").copy(
    #     Handle(BASE + "/pub/taxonomy/taxdump.tar.gz") \
    #     | "tar Oxf - names.dmp | grep scientific | cut -f1,3")

    print("Loading gene ...")
    db.table("gene").truncate()
    #db.table("gene").copy(
    Handle(BASE + "/gene/DATA/gene_info.gz") \
        | "sed 1d | grep -v '^#' | cut -f1-3,12 | awk '$1 != 366646'" \
        | "psql -d testdb -c \"COPY gene FROM STDIN NULL '-' \""
