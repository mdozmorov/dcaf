-- SELECT term.id, 
--        ARRAY(SELECT article.id FROM article WHERE article.abstract LIKE '%' || term.name || '%' LIMIT 5 )::int8[] as abstract
-- FROM (SELECT * FROM term LIMIT 5) as term;

SELECT DISTINCT term.id, ARRAY(article.id)
FROM (SELECT * FROM term LIMIT 10) as term
CROSS JOIN (SELECT * FROM article LIMIT 1000) AS article
WHERE article.abstract LIKE ('%' || term.name || '%')
GROUP BY term.id
ORDER BY term.id, article.id DESC
LIMIT 10;
