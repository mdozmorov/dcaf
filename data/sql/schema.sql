--DROP TABLE IF EXISTS taxon;
--DROP TABLE IF EXISTS assembly;
--DROP TABLE IF EXISTS contig;

--DROP TABLE IF EXISTS gene;
--DROP TABLE IF EXISTS feature;
--DROP TABLE IF EXISTS expression;
--DROP TABLE IF EXISTS sample;
--DROP TABLE IF EXISTS sample_expression;

--DROP TABLE IF EXISTS ontology;
--DROP TABLE IF EXISTS term;
--DROP TABLE IF EXISTS synonym;

--DROP TABLE IF EXISTS sample_term;
--DROP TABLE IF EXISTS feature_term;

CREATE TABLE taxon (
    id SERIAL PRIMARY KEY,
    name VARCHAR(255)
);

-- A particular genome assembly; e.g., GrCH37/hg19
CREATE TABLE IF NOT EXISTS assembly (
    id SERIAL PRIMARY KEY,
    taxon_id INTEGER,
    name VARCHAR(255),
    FOREIGN KEY (taxon_id) REFERENCES taxon (id),
    UNIQUE (name)
);

-- A chromosome or other contig belonging to a particular genome assembly.
CREATE TABLE IF NOT EXISTS contig (
    id SERIAL PRIMARY KEY,
    assembly_id INTEGER,
    name VARCHAR(255),
    size INTEGER,
    FOREIGN KEY (assembly_id) REFERENCES assembly(id)
);

CREATE TABLE IF NOT EXISTS gene (
    taxon_id INTEGER NOT NULL,
    id SERIAL PRIMARY KEY,
    symbol VARCHAR(255),
    name VARCHAR(1000),
    FOREIGN KEY (taxon_id) REFERENCES taxon(id)
);
CREATE INDEX gene_taxon_id_index ON gene(taxon_id);

CREATE TABLE IF NOT EXISTS feature (
    contig_id INTEGER, 
    span int8range
);

--------------------------------------
-- Expression dataset and sample types
--------------------------------------

CREATE TYPE continuous_imputable AS (
       value FLOAT,
       known BOOLEAN );
       
CREATE TYPE gender AS ENUM ('M', 'F');

CREATE TABLE IF NOT EXISTS platform (
    id SERIAL PRIMARY KEY,
    taxon_id INTEGER,
    accession VARCHAR(255),
    UNIQUE (accession)
);

CREATE TABLE IF NOT EXISTS sample (
    id SERIAL PRIMARY KEY,
    platform_id INTEGER,
    accession VARCHAR(255),
    title VARCHAR,
    characteristics VARCHAR,
    description VARCHAR,
    molecule VARCHAR,
    source VARCHAR,
    channel_count SMALLINT,

    CHECK(channel_count BETWEEN 1 AND 2),
    UNIQUE (accession),
    FOREIGN KEY (platform_id) REFERENCES platform (id)
);
CREATE INDEX ON sample (platform_id);

CREATE TABLE IF NOT EXISTS probe (
    id SERIAL PRIMARY KEY,
    platform_id INTEGER,
    accession VARCHAR(255),
    FOREIGN KEY (platform_id) REFERENCES platform (id)
);
CREATE INDEX ON probe(platform_id);

CREATE TABLE IF NOT EXISTS probe_gene (
    probe_id INTEGER,
    gene_id INTEGER,
    PRIMARY KEY (probe_id, gene_id),
    FOREIGN KEY (probe_id) REFERENCES probe (id),
    FOREIGN KEY (gene_id) REFERENCES gene (id)
);
CREATE INDEX ON probe_gene(gene_id);

CREATE TABLE IF NOT EXISTS expression (
    sample_id INTEGER,
    data float8[],
    known boolean[],
    FOREIGN KEY (sample_id) REFERENCES sample (id)
);
CREATE INDEX ON expression(sample_id);

CREATE OR REPLACE FUNCTION age(sample) RETURNS double precision AS $$
    SELECT cast((SELECT regexp_matches($1.characteristics, 
        'age: ([0-9]+\.*[0-9]*)'))[1] AS double precision)
$$ LANGUAGE SQL;

CREATE OR REPLACE FUNCTION gender(sample) RETURNS gender AS $$
    SELECT cast(upper((SELECT regexp_matches($1.characteristics, 
        '(?:gender|sex): ([MmFf])', 'i'))[1]) AS gender)
$$ LANGUAGE SQL;

-------------------------------
-- Ontology terms and traversal
-------------------------------

-- An OBO ontology, like GO, BRENDA, or SORD
CREATE TABLE IF NOT EXISTS ontology (
    id SERIAL PRIMARY KEY,
    namespace VARCHAR(255) NOT NULL,
    description VARCHAR(255) NOT NULL,
    UNIQUE (namespace)
);
INSERT INTO ontology VALUES (0, 'CORE', 'Core relations inherent to OBO format.');

-- An ontology term, denoting a concept like a disease,
--   phenotype, cellular component, or genomic element
CREATE TABLE IF NOT EXISTS term (
    id SERIAL PRIMARY KEY, 
    ontology_id INTEGER NOT NULL,
    accession VARCHAR(255),
    name VARCHAR(2000),
    UNIQUE (accession),
    FOREIGN KEY (ontology_id) REFERENCES ontology(id)
);
INSERT INTO term (id, ontology_id, name) VALUES (0, 0, 'is_a');
INSERT INTO term (id, ontology_id, name) VALUES (1, 0, 'part_of');

-- A synonym for a term. Each term may have many synonyms.
CREATE TABLE IF NOT EXISTS synonym (
    id SERIAL PRIMARY KEY,
    term_id INTEGER NOT NULL,
    synonym VARCHAR(2000),
    FOREIGN KEY (term_id) REFERENCES term(id)
);

-- Stores binary, directed relationships between terms
CREATE TABLE IF NOT EXISTS term_relation (
    agent INTEGER NOT NULL,
    target INTEGER NOT NULL,
    relation INTEGER NOT NULL,
    probability FLOAT DEFAULT 1,
    CHECK (probability BETWEEN 0 AND 1),
    FOREIGN KEY (agent) REFERENCES term (id),
    FOREIGN KEY (target) REFERENCES term(id),
    FOREIGN KEY (relation) REFERENCES term(id)
);

CREATE INDEX ON term_relation(agent);
CREATE INDEX ON term_relation(target);

--------------------------------
-- Sample and feature attributes
--------------------------------

CREATE TABLE IF NOT EXISTS sample_term (
    sample_id INTEGER NOT NULL,
    term_id INTEGER NOT NULL,
    probability FLOAT NOT NULL,
    UNIQUE (sample_id, term_id),
    FOREIGN KEY (sample_id) REFERENCES sample (id),
    FOREIGN KEY (term_id) REFERENCES term (id)
);
CREATE INDEX ON sample_term(sample_id);
CREATE INDEX ON sample_term(term_id);

CREATE TABLE IF NOT EXISTS gene_term (
    gene_id INTEGER NOT NULL,
    term_id INTEGER NOT NULL,
    probability FLOAT NOT NULL,
    UNIQUE (gene_id, term_id),
    FOREIGN KEY (gene_id) REFERENCES gene (id),
    FOREIGN KEY (term_id) REFERENCES term (id)
);
CREATE INDEX ON gene_term(gene_id);
CREATE INDEX ON gene_term(term_id);

--------------
-- Text mining
--------------

CREATE TABLE journal (
    id SERIAL PRIMARY KEY,
    issn VARCHAR,
    title VARCHAR
);

CREATE TABLE article (
    id SERIAL PRIMARY KEY,
    journal_id INTEGER,
    publication_date DATE,
    title VARCHAR,
    abstract VARCHAR,
    full_text VARCHAR,
    terms tsvector,
    FOREIGN KEY (journal_id) REFERENCES journal(id)
);

--CREATE TRIGGER article_tsvector_update BEFORE INSERT OR UPDATE
--       ON article FOR EACH ROW EXECUTE PROCEDURE
--       tsvector_update_trigger(terms, 'pg_catalog.english', title, abstract);
-- CREATE INDEX article_ix ON article USING gin(terms);

--CREATE TEXT SEARCH DICTIONARY iridescent (
--       TEMPLATE = thesaurus,
--       DictFile = iridescent
--);
