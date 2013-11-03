SET storage_engine=Aria;
SET NAMES 'utf8' COLLATE 'utf8_general_ci';

DROP TABLE IF EXISTS article;
CREATE TABLE article (
    id INTEGER NOT NULL,
    title TEXT,
    abstract TEXT,
    FULLTEXT(title),
    FULLTEXT(abstract),
    PRIMARY KEY(id)
);
