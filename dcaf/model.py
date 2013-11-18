"""
An API for SQL databases containing various genomic data.

This module contains utilities to build and query a database
containing:

* Expression profiles from GEO 
* Text articles from PubMed/MEDLINE
* Metadata from Entrez Gene and RefSeq

.. todo::
    * Add expression profiles from SRA
    * Various genomic feature tables from UCSC
    * Genomic coordinates for genes
"""

import logging
import os

from sqlalchemy import create_engine, Column, Integer, \
    String, ForeignKey, Float, Date
from sqlalchemy.orm import relationship
from sqlalchemy.schema import CheckConstraint, UniqueConstraint
from sqlalchemy.ext.declarative import as_declarative, declared_attr
from sqlalchemy.dialects.postgresql import INT8RANGE, ARRAY

from dcaf.io.medline import MedlineXMLFile

logging.basicConfig(level=logging.DEBUG,
                    format="%(asctime)s\t%(name)-12s\t%(levelname)-8s\t%(message)s",
                    datefmt="%m-%d %H:%M:%S")
log = logging.getLogger("dcaf.db")

def FK(column, index=True, nullable=False):
    return Column(Integer, ForeignKey(column, ondelete="CASCADE"),
                  index=index, nullable=nullable)

@as_declarative()
class Base(object):
    @declared_attr
    def __tablename__(cls):
        return "".join([x if x.islower() else "_"+x.lower()
                        for x in cls.__name__])[1:]

    id = Column(Integer, primary_key=True, nullable=False)

# Basic genome information

class Taxon(Base):
    name = Column(String)

    assemblies = relationship("Assembly", backref="taxon")
    genes = relationship("Gene", backref="taxon")
    
class Assembly(Base):
    taxon_id = FK("taxon.id")
    name = Column(String)

    contigs = relationship("Contig", backref="assembly")

class Contig(Base):
    assembly_id = FK("assembly.id")
    name = Column(String)
    size = Integer
    
    features = relationship("Feature", backref="contig")

class Gene(Base):
    taxon_id = FK("taxon.id")
    symbol = Column(String)
    name = Column(String)

class Feature(Base):
    contig_id = FK("contig.id")
    span = Column(INT8RANGE)

# Expression dataset & sample types

class Platform(Base):
    taxon_id = FK("taxon.id")
    accession = Column(String, unique=True)
    
    samples = relationship("Sample", backref="platform")

class Sample(Base):
    platform_id = FK("platform.id")
    accession = Column(String, unique=True)
    title = Column(String)
    characteristics = Column(String)
    description = Column(String)
    molecule = Column(String)
    source = Column(String)
    channel_count = Column(Integer)
    data = Column(ARRAY(Float, dimensions=1))
    
    _check = CheckConstraint("channel_count BETWEEN 1 AND 2")

# Ontology terms and traversal

class Ontology(Base):
    namespace = Column(String, unique=True)
    description = Column(String)

    terms = relationship("Term", cascade="all, delete, delete-orphan", 
                         backref="ontology")

class Term(Base):
    name = Column(String)
    accession = Column(String, unique=True)
    
    ontology_id = FK("ontology.id")
    synonyms = relationship("Synonym", cascade="all, delete, delete-orphan", 
                            backref="term")
    
    def __repr__(self):
        return "<Term %s: '%s'>" % (self.id, self.name)

class Synonym(Base):
    term_id = FK("term.id")
    synonym = Column(String)
    
class TermRelation(Base):
    agent = FK("term.id")
    target = FK("term.id")
    relation = FK("term.id")
    probability = Column(Integer)

    _check = CheckConstraint("probability BETWEEN 0 AND 1")
    _unique = UniqueConstraint("agent", "target", "relation")

# Sample and feature attributes

class SampleTerm(Base):
    sample_id = FK("sample.id")
    term_id = FK("term.id")
    probability = Column(Float)

    _check = CheckConstraint("probability BETWEEN 0 AND 1")
    _unique = UniqueConstraint("sample_id", "term_id")
    
class GeneTerm(Base):
    gene_id = FK("gene.id")
    term_id = FK("term.id")
    probability = Column(Float)

# Text mining

class Journal(Base):
    issn = Column(String)
    title = Column(String)
    
    articles = relationship("Article", backref="journal")

class Article(Base):
    journal_id = FK("journal.id")

    publication_date = Column(Date)
    title = Column(String)
    abstract = Column(String)
    full_text = Column(String)

###########################
# Database import functions
###########################

NCBI_BASE = "ftp://ftp.ncbi.nlm.nih.gov"

import dcaf.io
from plumbum.cmd import tar, grep, cut, zcat, awk, sed, head

# Use autocommit?

def initialize_db(session):
    log.info("Initializing DCAF database...")
    raw_db = session.connection().engine.raw_connection().connection

    # Add some built-in terms
    log.info("Adding core ontology terms...")
    session.add(Ontology(id=0, namespace="CORE",
                         description="Core relations inherent to OBO format."))
    session.add(Term(id=0, ontology_id=0, name="is_a"))
    session.add(Term(id=1, ontology_id=0, name="part_of"))
    session.commit()

    # Load NCBI taxa into 'taxon' table
    log.info("Loading NCBI taxa ...")
    url = NCBI_BASE + "/pub/taxonomy/taxdump.tar.gz"
    path = dcaf.io.download(url, return_path=True)

    c = raw_db.cursor()
    c.copy_from((tar["Oxzf", path, "names.dmp"] | grep["scientific"] \
                 | cut["-f1,3"]).popen().stdout, "taxon")
    raw_db.commit()

    # Load Entrez Gene IDs, names, symbols, etc, into 'gene'
    log.info("Loading NCBI Entrez Gene data ...")
    path = dcaf.io.download(NCBI_BASE + "/gene/DATA/gene_info.gz", 
                            return_path=True)
    c = raw_db.cursor()
    c.copy_from((zcat[path] | grep["-v", "^#"] | sed["1d"] \
                 | awk['$1 != 6617'] \
                 | awk['$1 != 543399'] \
                 | awk['BEGIN { OFS="\t" } $1 != 366646 { print $2,$1,$3,$12 }'] \
             ).popen().stdout, "gene")
    raw_db.commit()

def import_obo(session, handle, namespace, description):
    session.query(Ontology).filter(Ontology.namespace==namespace).delete()
    ontology = Ontology(namespace=namespace, description=description)
    session.add(ontology)

    is_a = session.query(Term).filter_by(name="is_a").first()
    part_of = session.query(Term).filter_by(name="part_of").first()
    relations = []

    obo = dcaf.io.OBOFile(handle)
    for term in obo:
        ontology.terms.append(Term(accession=term.id, name=term.name))
        for parent in term.is_a:
            relations.append((term.id, parent, is_a))
        for parent in term.part_of:
            relations.append((term.id, parent, part_of))
            
    for child_acc, parent_acc, rel in relations:
        child = session.query(Term).filter_by(accession=child_acc).first()
        parent = session.query(Term).filter_by(accession=parent_acc).first()
        session.add(TermRelation(agent=child.id, target=parent.id, 
                                 relation=rel.id, probability=1))
    session.commit()

def import_soft(session, path):
    parser = dcaf.io.SOFTParser(path)
    # FIXME: Make a better API

    taxon_id = parser._taxon_id
    platform_accession = parser._platform_accession
    
    platform = Platform(taxon_id=taxon_id, accession=platform_accession)
    session.add(platform)

    import itertools
    import numpy
    q = session.query(Gene).filter_by(taxon_id=taxon_id).with_entities(Gene.id)
    genes = sorted(itertools.chain(*q))

    for i,sample in enumerate(parser):
        v = sample["expression"]
        data = [v.get(g, numpy.nan) for g in genes]
        del sample["expression"]
        sample["data"] = data
        sample["platform_id"] = platform.id
        session.add(Sample(**sample))
        if (i % 100) == 0:
            session.commit()
            print("*", i, "records loaded.")
        session.commit()

def import_medline(session, folder):
    journal_ids = set()

    for file in os.listdir(folder):
        if file.endswith(".xml.gz"):
            path = os.path.join(folder, file)
            log.info("Importing MEDLINE articles from: %s" % path)
            with MedlineXMLFile(path) as handle:
                for article in handle:
                    journal = article.journal
                    if journal.id not in journal_ids:
                        session.add(Journal(
                            id=journal.id,
                            issn=journal.issn, 
                            title=journal.name))
                        journal_ids.add(journal.id)

                    session.add(Article(
                        journal_id=journal.id,
                        id=article.id,
                        publication_date=article.publication_date,
                        title=article.title,
                        abstract=article.abstract))
    session.commit()
                                
                    
if __name__ == "__main__":
    engine = create_engine("postgresql+psycopg2://wrendb/dcaf")
    db = engine.connect()
    Base.metadata.create_all(engine)

    from sqlalchemy.orm import sessionmaker
    Session = sessionmaker(bind=engine)

    session = Session()
    #import_medline(session, "/home/gilesc/data/text/")
    #initialize_db(session)

    #with open(dcaf.io.data("brenda.obo")) as h:
    #    import_obo(session, h, "BTO", "Brenda Tissue Ontology")

    #import_soft(session, "/home/gilesc/data/GEO/GPL96_family.soft.gz")
