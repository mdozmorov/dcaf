import logging

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

import dcaf.util

from .model import Base

__all__ = ["get_session"]

log = logging.getLogger("dcaf.db")

def get_session(**kwargs):
    # TODO: ignores user/password
    # TODO: memoize Session object
    cfg = dcaf.util.find_configuration()
    params = dict(list(cfg["database"].items()) + list(kwargs.items()))

    cstr = "postgresql+psycopg2://%(host)s:%(port)s/%(schema)s" % params
    engine = create_engine(cstr)
    db = engine.connect()
    Base.metadata.create_all(engine)
    Session = sessionmaker(bind=engine)
    return Session()
