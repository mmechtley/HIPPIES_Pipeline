from numpy import genfromtxt as _genfromtxt


def read_catalog(inName, keepCols=None, asrecord=True):
    """
    Read in a sextractor catalog and return as a dictionary
    """
    colNames = []
    ## Read in columns names from file header
    f = open(inName)
    for line in f:
        ## Stop when end of comments reached
        if not line.startswith('#'):
            break
        colNames += [line.split()[2]]
    f.close()

    if keepCols is not None:
        inds, colNames = zip(*[(i, name) for i, name in
                               enumerate(colNames) if name in keepCols])
    else:
        inds = range(len(colNames))

    ## Read in file data
    cat = _genfromtxt(inName, dtype=None, names=colNames, usecols=inds)
    if asrecord:
        return cat
    else:
        return cat.view(float).reshape(-1, len(colNames))
