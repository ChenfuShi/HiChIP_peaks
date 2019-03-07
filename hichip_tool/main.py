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
    ## here all inputs.


    folder = os.path.abspath("./../domain_caller/testdata/NaiveT_27ac_B1_T1")
    resfrag = os.path.abspath("./../domain_caller/testdata/MboI_resfrag_hg38.bed")
    sizes = None #or hicpro chromosome sizes
    temporary_loc = os.path.abspath("./../domain_caller/testdata")
    output_dir = os.path.abspath("./testdata")
    keeptemp = False
    threads=4
    FDR=0.01
    os.environ["OMP_NUM_THREADS"] = "threads" 
    os.environ["OPENBLAS_NUM_THREADS"] = "threads"
    os.environ["MKL_NUM_THREADS"] = "threads" 
    os.environ["VECLIB_MAXIMUM_THREADS"] = "threads" 
    os.environ["NUMEXPR_NUM_THREADS"] = "threads" 

    print("Info: \n HiC-Pro data folder: {} \n Restriction fragment file: {} \n Chromosome annotation file: {} \n Temporary location: {}".format(folder,resfrag,sizes,temporary_loc))
    print(" FDR: {}".format(FDR))
    print(" Output directory: {} \n Keep temporary files?: {} \n Threads(minimum is 4): {}".format(output_dir,keeptemp,threads))
    

    #apparently moving this should make sure that number of threads is respected in numpy?
    try:
        #this works only when installed
        from hichip_tool import interaction_to_sparse,sparse_to_peaks
    except:
        import interaction_to_sparse 
        import sparse_to_peaks
    
    CSR_mat,frag_index,frag_prop,frag_amount,valid_chroms,chroms_offsets = interaction_to_sparse.HiCpro_to_sparse(folder,resfrag,sizes,temporary_loc,keeptemp=keeptemp)


    smoothed_diagonal , refined_peaks = sparse_to_peaks.sparse_to_peaks(CSR_mat,frag_index,frag_prop,frag_amount,valid_chroms,chroms_offsets,output_dir,FDR=FDR,threads=threads,keeptemp=keeptemp)













if __name__=="__main__":
    main()