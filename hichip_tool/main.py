#########################################
# Author: Chenfu Shi
# Email: chenfu.shi@postgrad.manchester.ac.uk


# Main executor of the pipeline to call peaks from hicpro data

#########################################


def main():
    import os
    import argparse
    import logging
    import datetime
    import pickle
    ## here all inputs.


    parser = argparse.ArgumentParser(description="Peak calling from HiChIP data")

    parser.add_argument("-i", "--input", dest="hicpro_results",action="store",required=True,
                        help="HiC-Pro results directory containing validPairs file and others")
    parser.add_argument("-o", "--output", dest="output_directory",action="store",required=True,
                        help="Output directory")
    parser.add_argument("-r", "--resfrag", dest="resfrag",action="store",required=True,
                        help="HiCpro resfrag file")
    parser.add_argument("-p", "--prefix", dest="prefix",action="store",required=False, default=None,
                        help="Output file name prefix, if not provided will be name of HiC-Pro results directory")
    parser.add_argument("-f", "--FDR", dest="FDR",action="store",required=False, default=0.10, type=float,
                        help="False discovery rate, default = 0.10")                        
    parser.add_argument("-a", "--annotation", dest="sizes",action="store",required=False, default=None,
                        help="HiCpro chromosome annotation file, default uses human chromosomes, excludes chrY")
    parser.add_argument("-t", "--temporary_loc", dest="temporary_loc",action="store",required=False, default=None,
                        help="Temporary directory. If not supplied will be output directory")
    parser.add_argument("-w", "--worker_threads", dest="threads",action="store",required=False, default=4, type=int,
                        help="Number of threads, minimum 4")
    parser.add_argument("-k", "--keep_temp", dest="keeptemp",action="store_true", default=False,
                        help="Keep temporary files")
    parser.add_argument("-d", "--keep_diff", dest="keepdiff",action="store_true", default=False,
                        help="Prepare files for differential analysis")
    parser.add_argument("-s", "--offdiag", dest="off_diag",action="store",required=False, type=int, default=2,
                        help="How many off diagonal needs to be included")
    args = parser.parse_args()


    


    hicpro_results = os.path.abspath(args.hicpro_results)
    resfrag = os.path.abspath(args.resfrag)
    prefix=args.prefix
    sizes = args.sizes
    if args.temporary_loc == None:
        temporary_loc = os.path.abspath(args.output_directory)
    else:
        temporary_loc = os.path.abspath(args.temporary_loc)
    output_dir = os.path.abspath(args.output_directory)
    keeptemp = args.keeptemp
    keepdiff = args.keepdiff
    threads=args.threads
    FDR=args.FDR
    off_diag = args.off_diag
    if prefix == None:
        prefix = os.path.basename(hicpro_results)
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    logging.basicConfig(
        level=logging.INFO,
        format="%(levelname)s - %(message)s",
        handlers=[
        logging.FileHandler("{0}/{1}.log".format(output_dir, prefix + "log"), mode="a"),
        logging.StreamHandler()
    ]
    )    
    
    if threads < 4:
        threads = 4 
        logging.warning("Minimum threads is 4 !!")



    os.environ["OMP_NUM_THREADS"] = str(threads)
    os.environ["OPENBLAS_NUM_THREADS"] = str(threads)
    os.environ["MKL_NUM_THREADS"] = str(threads)
    os.environ["VECLIB_MAXIMUM_THREADS"] = str(threads)
    os.environ["NUMEXPR_NUM_THREADS"] = str(threads)

    logging.info("Welcome to my software!")
    logging.info(datetime.datetime.now())
    logging.info("Input variables")
    logging.info("HiC-Pro data folder: {} ".format(hicpro_results))
    logging.info("Restriction fragment file: {} ".format(resfrag))
    logging.info("Chromosome annotation file: {} ".format(sizes))
    logging.info("Temporary location: {} ".format(temporary_loc))
    logging.info("FDR: {} ".format(FDR))
    logging.info("Distance from diagonal included: {} ".format(off_diag))
    logging.info("Output directory: {} ".format(output_dir))
    logging.info("Output name prefix: {} ".format(prefix))
    logging.info("Keep temporary files?: {} ".format(keeptemp))
    logging.info("Threads(minimum is 4): {} ".format(threads))


    #apparently moving this should make sure that number of threads is respected in numpy?
    try:
        #this works only when installed
        from hichip_tool import interaction_to_sparse,sparse_to_peaks, quality_report
    except:
        import interaction_to_sparse 
        import sparse_to_peaks
        import quality_report
    
    CSR_mat,frag_index,frag_prop,frag_amount,valid_chroms,chroms_offsets = interaction_to_sparse.HiCpro_to_sparse(hicpro_results,resfrag,sizes,temporary_loc,prefix,keeptemp=keeptemp)


    smoothed_diagonal, refined_peaks ,quick_peaks, peak_p_vals , peaks_q_vals ,expected_background= sparse_to_peaks.sparse_to_peaks(CSR_mat,frag_index,frag_prop,frag_amount,valid_chroms,chroms_offsets,output_dir,prefix,off_diag,FDR=FDR,threads=threads,keeptemp=keeptemp)

    quality_report.quality_report(peak_p_vals,refined_peaks, smoothed_diagonal, output_dir, prefix)


    #if add an option to keep the data for the differential peak calling. then extra script that actually prepares the data for differential peak calling and goes into R
    #would still require the person to manually set design experiments and stuff.
    #save files for differential peak analysis
    if keepdiff == True:
        with open(os.path.join(output_dir,prefix + "diffpeak_data.pickle"),"wb") as picklefile:
            pickle.dump([smoothed_diagonal,refined_peaks,expected_background],picklefile)
    


    # # only for test and development purposes
    # with open(os.path.join(output_dir,prefix + "alldata.pickle"),"wb") as picklefile:
    #     pickle.dump([CSR_mat,frag_index,frag_prop,frag_amount,valid_chroms,chroms_offsets,smoothed_diagonal, refined_peaks ,quick_peaks, peak_p_vals , peaks_q_vals,expected_background],picklefile)




    logging.info(datetime.datetime.now())
    logging.info("Done!")

if __name__=="__main__":
    main()