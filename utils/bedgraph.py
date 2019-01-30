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

parser.add_argument("-i", "--input", dest="hicpro_results",action="store",required=True,
                    help="HiC-Pro results directory")
parser.add_argument("-o", "--output", dest="output_file",action="store",required=True,
                    help="Output file, will write temporary files in that directory")
parser.add_argument("-s", "--smoothing", dest="smoothing",action="store",required=False, type=int, Default=5,
                    help="Smoothing factor")
parser.add_argument("-d", "--offdiag", dest="off_diag",action="store",required=False, type=int, Default=2,
                    help="How many off diagonal needs to be included")
parser.add_argument("-r", "--resfrag", dest="resfrag",action="store",required=True,
                    help="HiCpro resfrag file")
parser.add_argument("-a", "--annotation", dest="valid_chroms",action="store",required=True,
                    help="Chromosomes annotation file")

args = parser.parse_args()

output_file = os.path.abspath(args.output_file)
temporary_loc = os.path.dirname(output_file)



CSR_mat,frag_index,frag_prop,frag_amount,valid_chroms,chroms_offsets,distribution_nice_fragments = HiCpro_to_sparse(args.hicpro_results,args.resfrag,args.valid_chroms,temporary_loc)

diagonal = extract_diagonal(CSR_mat,args.off_diag)
diagonal = moving_average(diagonal,args.smoothing)


# combine that with the fragment prop and you will have it :)