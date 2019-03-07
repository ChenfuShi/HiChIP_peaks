import gzip



def open_by_suffix(filename):
    if filename.endswith('.gz'):
        return gzip.open(filename, 'rt')
    else:
        return open(filename, 'r')