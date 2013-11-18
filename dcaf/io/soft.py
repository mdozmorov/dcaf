import collections
import gzip
import io
import itertools
import re

import pandas
import numpy

def as_float(item):
    try:
        return float(item)
    except:
        return numpy.nan
        
def read_table(lines, end_tag):
    buffer = io.StringIO()
    for line in lines:
        if line.startswith(end_tag):
            buffer.seek(0)
            return pandas.io.parsers.read_csv(buffer, sep="\t")
        buffer.write(line)

class Sample(object):
    """
    Represents a GEO GSM.
    """
    def __init__(self, id, expression, characteristics=[]):
        self.id = id
        self.expression = expression
        self.characteristics = characteristics

    def __repr__(self):
        return "<Sample %s with %d probes>" % (self.id, len(self.expression))
        
class SOFTParser(object):
    """
    A class to read from whole-platform SOFT datasets, such
    as those found at ftp://ftp.ncbi.nlm.nih.gov/geo/platforms/ . 

    It assumes the platform data will precede the samples.
    """
    def __init__(self, path):
        self._handle = gzip.open(path)
        self._lines = self._read_lines()

        while True:
            line = self._lines.__next__().strip()
            if line.startswith("!platform_table_begin"):
                break
            elif line.startswith("!Platform_taxid"):
                self._taxon_id = int(self._read_value(line))
            elif line.startswith("!Platform_geo_accession"):
                self._platform_accession = self._read_value(line)

        platform = self._read_table("!platform_table_end")

        if "Entrez_Gene_ID" in platform.columns:
            columns = list(platform.columns)
            columns[columns.index("Entrez_Gene_ID")] = "ENTREZ_GENE_ID"
            platform.columns = columns

        if "GENE" in platform.columns:
            columns = list(platform.columns)
            columns[columns.index("GENE")] = "ENTREZ_GENE_ID"
            platform.columns = columns

        if not "ENTREZ_GENE_ID" in platform.columns: 
            raise Exception("Cannot handle this platform as there is \
            no mapping of probes to Entrez Gene IDs.")

        self._probe_to_gene = {}
        for probe, gene in zip(platform["ID"], platform["ENTREZ_GENE_ID"]):
            # Take only probes that map to exactly one gene
            # Those that don't will have '///' in the name
            try:
                self._probe_to_gene[probe] = int(gene)
            except ValueError:
                continue

    def _read_lines(self):
        for line in self._handle:
            try:
                yield line.decode("UTF-8")
            except UnicodeDecodeError:
                continue
                
    def _read_table(self, end_tag):
        """Read a platform or series table into a DataFrame."""
        buffer = io.StringIO()
        for line in self._lines:
            if line.startswith(end_tag):
                buffer.seek(0)
                return pandas.io.parsers.read_csv(buffer, sep="\t")
            buffer.write(line)

    def _read_until(self, tag):
        """Advance the internal lines iterator to the given tag."""
        for line in self._lines:
            if line.startswith(tag):
                return

    def _read_value(self, line):
        return "".join(line.strip().split(" = ")[1:])

    def __iter__(self):
        """Iterate through the Sample objects in this SOFT file."""
        #sample = {"accession" : sample_accession}
        sample = {}
        for line in self._lines:
            if line.startswith("!Sample"):
                try:
                    key, value = line.strip().split(" = ", 1)
                    key = key.replace("!Sample_", "")
                    if key == "characteristics_ch1":
                        sample.setdefault("characteristics",[])
                        sample["characteristics"].append(value)
                    elif key == "title":
                        sample["title"] = value
                    elif key == "description":
                        sample.setdefault("description", [])
                        sample["description"].append(value)
                    elif key == "channel_count":
                        sample["channel_count"] = int(value)
                    elif key == "molecule_ch1":
                        sample["molecule"] = value
                    elif key == "source_name_ch1":
                        sample["source"] = value
                except Exception as e:
                    continue
            elif line.startswith("!sample_table_begin"):
                try:
                    expression = read_table(self._lines, "!sample_table_end")
                    expression_vector = pandas.Series(
                        [as_float(item) for item in expression["VALUE"]],
                        index=expression["ID_REF"]).groupby(self._probe_to_gene)\
                                              .max()
                    sample["expression"] = expression_vector
                    yield sample
                except Exception as e:
                    print(e)
                    continue
 
def parseSOFT(path):
   parser = SOFTParser(path) 
   return ExpressionSet(parser.phenotype_data,
                        parser.feature_data, parser.expression)
