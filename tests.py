from scripts.interaction_to_sparse import *
from scripts.sparse_to_peaks import *

import pickle
# things needed

# hicpro directory containing all the valid ints files
# restriction fragment definitions
# sizes of chromosomes just to get the right chromosomes

folder = os.path.abspath("/mnt/iusers01/jw01/mdefscs4/psa_functional_genomics/OurData_Hichip/data/homebrew-HiC-Pro/hic_results/data/HaCat_27ac_st")
resfrag = os.path.abspath("./../domain_caller/testdata/MboI_resfrag_hg38.bed")
sizes = os.path.abspath("./annotations/hg38.txt")
temporary_loc = os.path.abspath("./testdata/combined_Mumbach")

CSR_mat,frag_index,frag_prop,frag_amount,valid_chroms,chroms_offsets = HiCpro_to_sparse(folder,resfrag,sizes,temporary_loc)
smoothed_diagonal, refined_peaks = sparse_to_peaks(CSR_mat,frag_index,frag_prop,frag_amount,valid_chroms,chroms_offsets,temporary_loc,FDR=0.10)

scipy.sparse.save_npz('./testdata/sparse_matrix_HACAT_st.npz', CSR_mat)

# with open("./testdata/variables.pi","wb") as picklefile:
#    pickle.dump([frag_index,frag_prop,frag_amount,valid_chroms,chroms_offsets],picklefile)


    