#########################################
# Author: Chenfu Shi
# Email: chenfu.shi@postgrad.manchester.ac.uk


# Function inputs a HiChIP CSR_mat and extracts the peaks. 
# Extracts short range interactions from these peaks, forms a vector. Calls significantly enriched fragments then cleans up and merges peaks together. 
# Returns list of peaks as fragment numbers

#########################################

import scipy
import scipy.sparse
import numpy
import os
import re
import multiprocessing
import subprocess

def sparse_to_peaks(CSR_mat):
    """Wrapper function to call individual funcitons"""

























if __name__=="__main__":
    CSR_mat = scipy.sparse.load_npz('./testdata/sparse_matrix.npz')