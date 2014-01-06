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
