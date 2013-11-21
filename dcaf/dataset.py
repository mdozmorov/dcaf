from pandas import read_excel, Categorical, DataFrame

import dcaf.io

def _load_hfd_hippocampus():
    cols = [0]
    for i in range(4):
        cols.extend(range(i*5+1+i, i*5+6+i))

    path = dcaf.io.data("expression/hfd_hippocampus.xlsx")
    X = read_excel(path, "AHFD hippocampi", 
                   na_values="Undetermined", 
                   parse_cols=cols, skiprows=4, index_col=0)
    X.index = [x.split("-")[0] for x in X.index]
    #baseline = X.head(3).prod().pow(1/3)
    #X = ((2 ** baseline) / (2 ** X)).T
    X = X.T

    P = DataFrame(
        {"Age":
         Categorical((["Young"] * 10) + (["Old"] * 10)),
         "Diet":
         Categorical(((["Control"] * 5) + (["HFD"] * 5)) * 2)
        }, index=X.index)

    return X, P

class DatasetNotFound(Exception):
    def __init__(self, name):
        msg = "No loader for dataset named '%s'." % name
        super(DatsetNotFound, self).__init__(msg)

def load_dataset(name):
    if name == "hfd_hippocampus":
        return _load_hfd_hippocampus()
    raise DatasetNotFound(name)
