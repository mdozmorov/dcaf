# Notes:
# - http://www.postgresql.org/docs/9.1/static/intarray.html
# - load with 'CREATE EXTENSION intarray'

# TODO:
# - Load genes, DO, CHEBI
# - Handle synonyms correctly
# - Build web interface

from itertools import groupby, chain
from operator import itemgetter

from acora import AcoraBuilder

class TermSearch(object):
    def __init__(self):
        self._builder_cs = AcoraBuilder()
        self._builder_ci = AcoraBuilder()
        self._text_to_id = {}
    
    def add(self, id, term, case_sensitive=False):
        builder = self._builder_cs if case_sensitive else self._builder_ci
        if not case_sensitive:
            term = term.lower()
        builder.add(term)
        self._text_to_id[term] = id

    def build(self):
        self._cs = self._builder_cs.build()
        self._ci = self._builder_ci.build()

    def search(self, text):
        # FIXME: prioritize CS over CI before merge
        #matches = chain(self._cs.finditer(text), self._ci.finditer(text.lower()))
        matches = self._ci.finditer(text.lower())
        for pos, match_set in groupby(matches, itemgetter(1)):
            item, start = max(match_set)
            yield (self._text_to_id[item], start, start+len(item))

from collections import defaultdict
            
import dcaf
from dcaf.db.model import *

def find_postings_for_ontology_terms(session):
    ts = TermSearch()
    term_name = {}
    for term in session.query(Term).limit(1000):
        ts.add(term.id, term.name)
        for synonym in term.synonyms:
            ts.add(term.id, synonym.synonym)
        term_name[term.id] = term.name
    ts.build()
    print("Built term search automaton...")

    index = defaultdict(set)

    for i,article in enumerate(session.query(Article).limit(500000)):
        text = " ".join([article.title, article.abstract or ""]) 
        for term_id, start, end in ts.search(text):
            index[term_id].add(article.id)

def sync_postings():
    session = dcaf.db.get_session()
    postings = find_postings_for_ontology_terms(session)
    for term_id, articles in index.items():
        session.add(Posting(term_id=term_id, articles=articles))
    session.commit()

def find_postings_for_text(text):
    """
    FIXME: Search using tsvector
    FIXME: Search case insensitive (sqlalchemy.func.lower)
    """
    session = dcaf.db.get_session()
    text = text.lower()
    postings = session.query(Article.id)\
                      .filter((Article.title.contains(text)) | 
                              (Article.abstract.contains(text)))
    return set(r[0] for r in postings)

import pandas

def similar_terms(text):
    rows = []
    session = dcaf.db.get_session()
    postings = find_postings_for_text(text)
    for term_id, term_name, term_postings in \
        session.query(Term.id, Term.name, Posting.articles)\
               .join(Posting):
        term_postings = set(term_postings)
        jaccard = len(postings & term_postings) / \
                  len(postings | term_postings)
        if jaccard:
            rows.append((term_name, jaccard))
    return pandas.DataFrame.from_records(rows, 
                                         columns=["Term", "Jaccard Coefficient"])
        
            
from flask import Flask, request
from mako.template import Template
from mako.lookup import TemplateLookup

app = Flask(__name__)
lookup = TemplateLookup(directories=["template"])

def render_table(df):
    tmpl = lookup.get_template("table.html")
    header = []
    rows = []
    if df is not None:
        header = list(map(str, df.columns))
        rows = [list(row)[1:] for row in df.to_records()]
    return tmpl.render(rows=rows, header=header)

@app.route("/")
def main():
    term = request.args.get("term")
    df = None
    if term:
        df = similar_terms(term).sort("Jaccard Coefficient", 
                                      ascending=False)
    tmpl = Template(open("template/term_search.html").read())
    return tmpl.render(result=render_table(df))
    
if __name__ == "__main__":
    app.run()
