#########################################
# Author: Chenfu Shi
# Email: chenfu.shi@postgrad.manchester.ac.uk


# Main executor of the pipeline

#########################################


# do one

# check results
# print results such as how many reads

# do second

# repeat

def main():
    import os
    try:
        #this works only when installed
        from hichip_tool import interaction_to_sparse,sparse_to_peaks
    except:
        import interaction_to_sparse 
        import sparse_to_peaks
    ## here all inputs.


    folder = os.path.abspath("./../domain_caller/testdata/NaiveT_27ac_B1_T1")
    resfrag = os.path.abspath("./../domain_caller/testdata/MboI_resfrag_hg38.bed")
    sizes = None #or hicpro chromosome sizes
    temporary_loc = os.path.abspath("./../domain_caller/testdata")
    output_dir = os.path.abspath("./testdata")
    keeptemp = False




    CSR_mat,frag_index,frag_prop,frag_amount,valid_chroms,chroms_offsets = interaction_to_sparse.HiCpro_to_sparse(folder,resfrag,sizes,temporary_loc,keeptemp=keeptemp)


    smoothed_diagonal , refined_peaks = sparse_to_peaks.sparse_to_peaks(CSR_mat,frag_index,frag_prop,frag_amount,valid_chroms,chroms_offsets,output_dir)













if __name__=="__main__":
    main()