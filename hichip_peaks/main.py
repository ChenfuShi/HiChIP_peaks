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
    
    ## parse inputs from command line
    parser = argparse.ArgumentParser(description="Peak calling from HiChIP data")

    parser.add_argument("-i", "--input", dest="hicpro_results",action="store",required=True,
                        help="HiC-Pro results directory containing validPairs file and others")
    parser.add_argument("-o", "--output", dest="output_directory",action="store",required=True,
                        help="Output directory")
    parser.add_argument("-r", "--resfrag", dest="resfrag",action="store",required=True,
                        help="HiCpro resfrag file")
    parser.add_argument("-p", "--prefix", dest="prefix",action="store",required=False, default=None,
                        help="Output file name prefix, if not provided will be name of HiC-Pro results directory")
    parser.add_argument("-f", "--FDR", dest="FDR",action="store",required=False, default=0.01, type=float,
                        help="False discovery rate, default = 0.01")                        
    parser.add_argument("-a", "--annotation", dest="sizes",action="store",required=False, default=None,
                        help="HiCpro chromosome annotation file, default uses human chromosomes, excludes chrY")
    parser.add_argument("-t", "--temporary_loc", dest="temporary_loc",action="store",required=False, default=None,
                        help="Temporary directory. If not supplied will be output directory")
    parser.add_argument("-w", "--worker_threads", dest="threads",action="store",required=False, default=4, type=int,
                        help="Number of threads, minimum 4. Warning: Increasing this significantly increases RAM usage")
    parser.add_argument("-k", "--keep_temp", dest="keeptemp",action="store_true", default=False,
                        help="Keep temporary files")
    parser.add_argument("-d", "--keep_diff", dest="keepdiff",action="store_true", default=False,
                        help="Prepare files for differential analysis")
    parser.add_argument("-s", "--offdiag", dest="off_diag",action="store",required=False, type=int, default=2,
                        help="How many off diagonal needs to be included (default = 2)")
    parser.add_argument("-x", "--chromX", dest="chromX",action="store_true",required=False, default=False,
                        help="Want to compensate Sex chromosomes weights? Requires specify annotation(SIZES) containing chrX and chrY")
    parser.add_argument("-c", "--class_store", dest="class_store",action="store_true",required=False, default=False,
                        help="Store sparse site_matrix object for further use")
    args = parser.parse_args()

    # parse arguments and check sanity of inputs
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
    chromX = args.chromX
    class_store = args.class_store
    if prefix == None:
        prefix = os.path.basename(hicpro_results)
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    
    # Create logging file. prints both to screen and to log file
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

    # set number of threads for numpy
    os.environ["OMP_NUM_THREADS"] = str(threads)
    os.environ["OPENBLAS_NUM_THREADS"] = str(threads)
    os.environ["MKL_NUM_THREADS"] = str(threads)
    os.environ["VECLIB_MAXIMUM_THREADS"] = str(threads)
    os.environ["NUMEXPR_NUM_THREADS"] = str(threads)

    # import functions that do actual analysis
    # apparently moving this after the threads set should make sure that number of threads is respected in numpy?
    try:
        #this works only when installed, which is the large majority of use cases
        from hichip_peaks import interaction_to_sparse,sparse_to_peaks, quality_report
        from hichip_peaks.site_matrix_class import site_matrix
        from hichip_peaks.__init__ import __version__
    except:
        import interaction_to_sparse 
        import sparse_to_peaks
        import quality_report
        from __init__ import __version__

    # log inputs
    logging.info("Welcome to HiChIP-Peaks!")
    logging.info("Version: {}".format(__version__))
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
    logging.info("Sex chromosomes correction: {} ".format(chromX))
    logging.info("Store object containing sparse matrix: {} ".format(class_store))

    # Start of actualy analysis
    # call first function to create sparse matrix representation of hichip data

    CSR_mat,frag_index,frag_prop,frag_amount,valid_chroms,chroms_offsets = interaction_to_sparse.HiCpro_to_sparse(hicpro_results,resfrag,sizes,temporary_loc,prefix,keeptemp=keeptemp)

    # call peak calling algorithm. This also saves peaks and bedgraph

    smoothed_diagonal, refined_peaks ,quick_peaks, peak_p_vals , peaks_q_vals ,expected_background= sparse_to_peaks.sparse_to_peaks(CSR_mat,frag_index,frag_prop,frag_amount,valid_chroms,chroms_offsets,output_dir,prefix,off_diag,chromX,FDR=FDR,threads=threads,keeptemp=keeptemp)

    # call quality report generator
    logging.info("#######################################")
    logging.info("Creating report with summary statistics")
    quality_report.quality_report(peak_p_vals,refined_peaks, smoothed_diagonal, output_dir, prefix)

    # save files for differential peak analysis
    if keepdiff == True:
        logging.info("#######################################")
        logging.info("Saving file for differential peak analysis")
        with open(os.path.join(output_dir,prefix + "diffpeak_data.pickle"),"wb") as picklefile:
            pickle.dump([smoothed_diagonal,refined_peaks,expected_background],picklefile)
    
    if class_store == True:
        logging.info("#######################################")
        logging.info("Storing sparse matrix object for future analysis")
        mat_obj = site_matrix(CSR_mat, valid_chroms,frag_prop,frag_index,frag_amount,chroms_offsets)
        peaks = mat_obj.register_peaks(refined_peaks)
        with open(os.path.join(output_dir,prefix + "mat_obj.pickle"),"wb") as picklefile:
            pickle.dump(mat_obj,picklefile)

    # # only for test and development purposes
    # with open(os.path.join(output_dir,prefix + "alldata.pickle"),"wb") as picklefile:
    #     pickle.dump([CSR_mat,frag_index,frag_prop,frag_amount,valid_chroms,chroms_offsets,smoothed_diagonal, refined_peaks ,quick_peaks, peak_p_vals , peaks_q_vals,expected_background],picklefile)

    logging.info("#######################################")
    logging.info(datetime.datetime.now())
    logging.info("Done!")


if __name__=="__main__":
    main()