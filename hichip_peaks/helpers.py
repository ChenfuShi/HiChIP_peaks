#########################################
# Author: Chenfu Shi
# Email: chenfu.shi@postgrad.manchester.ac.uk


# some helper functions

#########################################


import gzip
import scipy
import scipy.sparse
import scipy.stats , scipy.interpolate
import statsmodels.stats.multitest , statsmodels , statsmodels.api
import math
import numpy
import os
import re
import multiprocessing
import subprocess
import matplotlib.pyplot
import itertools
import logging
import pickle


def open_by_suffix(filename):
    if filename.endswith('.gz'):
        return gzip.open(filename, 'rt')
    else:
        return open(filename, 'r')
















### functions for sparse matrix format





