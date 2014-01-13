from flask import Flask, request
from mako.template import Template
from mako.lookup import TemplateLookup

from dcaf.text import similar_terms

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

@app.route("/literature")
def literature_main():
    """
    Main page for literature analysis of terms.
    """
    term = request.args.get("term")
    df = None
    if term:
        df = similar_terms(term).sort("Jaccard Coefficient", 
                                      ascending=False)
    tmpl = lookup.get_template("term_search.html")
    return tmpl.render(result=render_table(df))
    
@app.route("/")
def main():
    tmpl = lookup.get_template("index.html")
    return tmpl.render()

if __name__ == "__main__":
    app.run()
