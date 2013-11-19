"""
Command-line interface for manipulating the ``dcaf`` PostgreSQL database
(mostly data import). 
"""

import argparse

import dcaf.util

from .load import import_medline, import_soft, import_obo, initialize_db
from .core import get_session, log

@dcaf.util.entry_point
def dcafdb(argv):
    """
    CLI entry point.
    """
    p = argparse.ArgumentParser(prog="dcafdb")

    # Parsers for subcommands
    sp = p.add_subparsers()
    
    sp_initialize_db = sp.add_parser("initialize-db")
    sp_initialize_db.set_defaults(func=lambda session,args: 
                                  initialize_db(session))

    sp_import_medline = sp.add_parser("import-medline")
    sp_import_medline.add_argument("path")
    sp_import_medline.set_defaults(func=lambda session,args: 
                                   import_medline(session, args.path))
    
    sp_import_soft = sp.add_parser("import-soft")
    sp_import_soft.add_argument("path")
    sp_import_soft.set_defaults(func=lambda session, args:
                                import_soft(session, args.path))
    
    sp_import_obo = sp.add_parser("import-obo")
    sp_import_obo.add_argument("path")
    sp_import_obo.add_argument("namespace")
    sp_import_obo.add_argument("description")
    sp_import_obo.set_defaults(func=lambda session, args:
                               import_obo(session, args.path, 
                                          args.namespace,
                                          args.description))

    args = p.parse_args(argv)

    # TODO: Allow user to specify DB connection parameters

    if "func" in args:
        session = get_session()
        from .load import _import_go
        args.func(session, args)

        
if __name__ == "__main__":
    dcafdb()
