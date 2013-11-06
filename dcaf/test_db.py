def test_connect():
    import dcaf.db

    db = dcaf.db.DCAFConnection.from_configuration()
