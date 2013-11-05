import dcaf.db

def test_connect():
    db = dcaf.db.Connection(database="gfa", host="wrendb")
