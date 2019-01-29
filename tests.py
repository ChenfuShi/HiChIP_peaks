from scripts.interaction_to_sparse import *

import pickle
# things needed

# hicpro directory containing all the valid ints files
# restriction fragment definitions
# sizes of chromosomes just to get the right chromosomes

folder = os.path.abspath("/mnt/iusers01/jw01/mdefscs4/psa_functional_genomics/combined_hichip/data/HiC-pro_results/hic_results/data/mumbach_27ac")
resfrag = os.path.abspath("./../domain_caller/testdata/MboI_resfrag_hg38.bed")
sizes = os.path.abspath("./annotations/hg38.txt")
temporary_loc = os.path.abspath("./testdata/combined_Mumbach")

CSR_mat,frag_index,frag_prop,frag_amount,valid_chroms,chroms_offsets,distribution_nice_fragments = HiCpro_to_sparse(folder,resfrag,sizes,temporary_loc)

scipy.sparse.save_npz('./testdata/sparse_matrix_mumbach.npz', CSR_mat)

with open("./testdata/variables.pi","wb") as picklefile:
    pickle.dump([frag_index,frag_prop,frag_amount,valid_chroms,chroms_offsets,distribution_nice_fragments],picklefile)


    