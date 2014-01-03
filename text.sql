-- SELECT term.id, 
--        ARRAY(SELECT article.id FROM article WHERE article.abstract LIKE '%' || term.name || '%' LIMIT 5 )::int8[] as abstract
-- FROM (SELECT * FROM term LIMIT 5) as term;

CREATE EXTENSION IF NOT EXISTS intarray;

CREATE MATERIALIZED VIEW term_profile (term_id, articles) AS
    SELECT term.id AS term_id, uniq(sort(array_agg(article.id))) as articles
    FROM term --(SELECT * FROM term LIMIT 100) as term
    CROSS JOIN (SELECT * FROM article LIMIT 10000) AS article
    WHERE 
          (article.abstract LIKE ('%' || term.name || '%') OR
                            article.title LIKE ('%' || term.name || '%'))
    GROUP BY term.id;
