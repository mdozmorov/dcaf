""" 
Synchronize the database with external data sources using the
configured settings. 
"""

from ..io import generic_open as open
from .model import *
from . import get_session

class DataSource(object):
    def download(self):
        """
        Download the raw data files for this DataSource to 
        on-disk cache.
        """
        pass
    
    def clear(self, session):
        """
        Delete all data from this DataSource in the database.
        """
        
    def is_loaded(self, session):
        """
        Check if this DataSource has already been loaded.
        """
        pass

    def sync(self, session):
        """
        Load or update the database with the latest version 
        of this DataSource.
        """
    
class OBODataSource(DataSource):
    """
    A DataSource backed by one of the Open Biomedical Ontologies.
    """
    def __init__(self, namespace, description, url):
        self._namespace = namespace
        self._description = description
        self._url = url
    
    def clear(self, session):
        session.query(Ontology).filter(Ontology.namespace==self._namespace).delete()
    
    def is_loaded(self, session):
        return session.query(Ontology)\
                      .filter(Ontology.namespace==self._namespace)\
                      .first() is not None
    
    def sync(self, session):
        if self.is_loaded:
            return

        
        with open(self._url) as handle:
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

def main():
    session = get_session()
    obo = OBODataSource("GO", "Gene Ontology",
                        "http://www.geneontology.org/ontology/obo_format_1_2/gene_ontology_ext.obo")
    obo.clear()
    print(obo.is_loaded())
    obo.sync(session)
    print(obo.is_loaded(session))

if __name__ == "__main__":
    main()
