"""
Query the dcaf database.

TODO:

- build an Ontology object
- 
"""

from .model import *

def fetch_expression(session, taxon_id=9606):
    q = session.query(Sample)\
               .join(Platform.samples)\
               .filter(Platform.taxon_id==taxon_id)\
               .filter(Sample.characteristics != None)\
               .limit(5)

    for item in q:
        print(item)
        #if item.age:
        #    print(item.accession, item.characteristics, item.age)
        #print(item.id)
        #print(item.accession)
        #print(item.title)
        #print(item.age)
        #print(item.description)

if __name__ == "__main__":
    from .core import get_session

    session = get_session()
    fetch_expression(session)
