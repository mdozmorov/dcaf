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

import os
import re

from sqlalchemy import create_engine, Column, Integer, \
    String, ForeignKey, Float, Date
from sqlalchemy.orm import relationship, deferred
from sqlalchemy.schema import CheckConstraint, UniqueConstraint
from sqlalchemy.ext.declarative import as_declarative, declared_attr
from sqlalchemy.ext.hybrid import hybrid_property

from sqlalchemy.dialects.postgresql import INT8RANGE, ARRAY, TEXT


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
    data = deferred(Column(ARRAY(Float, dimensions=1)))
    
    @property
    def age(self):
        ch = self.characteristics
        if ch:
            m = re.search('[Aa]ge[A-Za-z ()]*: ([0-9]+\.*[0-9]*)', ch)
            if m:
                return float(m.group(1))

    _check = CheckConstraint("channel_count BETWEEN 1 AND 2")

# Ontology terms and traversal

class Ontology(Base):
    namespace = Column(String, unique=True)
    description = Column(String)

    terms = relationship("Term", cascade="all, delete, delete-orphan", 
                         backref="ontology")

class Term(Base):
    name = Column(String)
    description = Column(String)
    long_description = Column(String)
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

#############
# Text mining
#############

from sqlalchemy import types
from sqlalchemy.ext.compiler import compiles

class tsvector(types.TypeDecorator):
    impl = types.UnicodeText

@compiles(tsvector, "postgresql")
def compile_tsvector(element, compiler, **kwargs):
    return "tsvector"

class Journal(Base):
    issn = Column(String)
    title = Column(String)
    
    articles = relationship("Article", backref="journal")
    
class Article(Base):
    journal_id = FK("journal.id")

    publication_date = Column(Date)
    title = Column(TEXT)
    abstract = Column(TEXT)
    full_text = Column(TEXT)

class Posting(Base):
    term_id = FK("term.id")
    
    articles = Column(ARRAY(Integer))
    #title = Column(ARRAY(Integer))
    #abstract = Column(ARRAY(Integer))
    #sentence = Column(ARRAY(Integer)) 
    #full_text = Column(ARRAY(INTEGER))
