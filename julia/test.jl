using DBI
using PostgreSQL

connect(Postgres, "wrendb", "gilesc", "", "dcaf", 5432) do db
    stmt = prepare(db, "SELECT id, data FROM sample LIMIT 5")
    X = zeros(Float64, (5,5))
    for (i,row) in enumerate (execute (stmt))
        X[i,:] = row[2][1:5]
    end
    println (X)
    finish (stmt)
end
