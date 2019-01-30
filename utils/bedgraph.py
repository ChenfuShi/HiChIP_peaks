#########################################
# Author: Chenfu Shi
# Email: chenfu.shi@postgrad.manchester.ac.uk


# makes bedgraphs from scratch


#########################################

import os, sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from scripts.interaction_to_sparse import *
from scripts.sparse_to_peaks import *

import argparse

parser = argparse.ArgumentParser(description="input HiCPro results and generates a bedgraph file")

parser.add_argument("-i")


HiCpro_to_sparse()

extract_diagonal()

# combine that with the fragment prop and you will have it :)