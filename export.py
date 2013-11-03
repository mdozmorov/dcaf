import grtk.db

db = grtk.db.connect()
c = db.cursor()
c.execute("SELECT id FROM gene WHERE taxon_id=9606 ORDER BY id")
print(*["accession"] + [row[0] for row in c], sep="\t")
c = db.cursor("export")
c.itersize = 1
c.execute("""
SELECT sample.accession,data FROM expression 
INNER JOIN sample
ON sample.id=expression.sample_id""")
for row in c:
    print(*[row[0]] + row[1], sep="\t")
