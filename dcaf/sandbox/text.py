# Notes:
# - http://www.postgresql.org/docs/9.1/static/intarray.html
# - load with 'CREATE EXTENSION intarray'

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

session = dcaf.db.get_session()

ts = TermSearch()
term_name = {}
for term in session.query(Term).limit(20):
    ts.add(term.id, term.name)
    term_name[term.id] = term.name
ts.build()

index = defaultdict(set)

for i,article in enumerate(session.query(Article).limit(500000)):
    text = " ".join([article.title, article.abstract or ""]) 
    for term_id, start, end in ts.search(text):
        index[term_id].add(article.id)

for term_id, articles in index.items():
    print(term_id, term_name[term_id], len(articles))
