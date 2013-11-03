CREATE OR REPLACE FUNCTION age(sample) RETURNS double precision AS $$
    SELECT cast((SELECT regexp_matches($1.characteristics, 
        '[Aa]ge[A-Za-z ()]+: ([0-9]+\.*[0-9]*)'))[1] AS double precision)

$$ LANGUAGE SQL;

CREATE OR REPLACE FUNCTION gender(sample) RETURNS gender AS $$
    SELECT cast(upper((SELECT regexp_matches($1.characteristics, 
        '(?:gender|sex): ([MmFf])', 'i'))[1]) AS gender)
$$ LANGUAGE SQL;

CREATE OR REPLACE FUNCTION ontology_links(namespace VARCHAR) RETURNS TABLE (child INTEGER, parent INTEGER) AS $$
    SELECT agent, target FROM term_relation r
    INNER JOIN term t1 ON r.agent=t1.id
    INNER JOIN term t2 ON r.target=t2.id
    INNER JOIN ontology o1 ON t1.ontology_id=o1.id
    INNER JOIN ontology o2 ON t2.ontology_id=o2.id
    WHERE r.relation=0
    AND o1.namespace=$1 AND o2.namespace=$1
$$ LANGUAGE SQL;

--CREATE OR REPLACE FUNCTION gene_expression(taxon_id INTEGER)
--    RETURNS TABLE(data float8[]) AS $$
    --SELECT --id as gene_id, 
--    SELECT array_agg(x) FROM 
--        (SELECT data[ix] AS x FROM 
--            (SELECT id, row_number() OVER (ORDER BY id) AS ix FROM 
--                gene WHERE taxon_id=$1) AS q) AS foo
--$$ LANGUAGE SQL;

DROP VIEW foo;
CREATE OR REPLACE VIEW foo AS
    SELECT id1, id2
        FROM
        (SELECT row_number() OVER (ORDER BY id) AS ix1, 
            id AS id1
            FROM gene WHERE taxon_id=9606) A
        CROSS JOIN
        (SELECT row_number() OVER (ORDER BY id) AS ix2, 
            id AS id2 FROM gene WHERE taxon_id=9606) B
    LIMIT 10
