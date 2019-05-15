#########################################
# Author: Chenfu Shi
# Email: chenfu.shi@postgrad.manchester.ac.uk


# makes bedgraphs from scratch
# bedgraphs are generated from the peak calling. also these bedgraphs are slightly different. you can set smoothing factor and offdiag factor separately.

#########################################
def main():
    import os, sys

    try:
        #this works only when installed
        from hichip_peaks import interaction_to_sparse,sparse_to_peaks
    except:
        import interaction_to_sparse 
        import sparse_to_peaks
    import argparse
    import math

    parser = argparse.ArgumentParser(description="input HiCPro results and generates a bedgraph file")

    parser.add_argument("-i", "--input", dest="hicpro_results",action="store",required=True,
                        help="HiC-Pro results directory")
    parser.add_argument("-o", "--output", dest="output_file",action="store",required=True,
                        help="Output file, will write temporary files in that directory")
    parser.add_argument("-s", "--smoothing", dest="smoothing",action="store",required=False, type=int, default=3,
                        help="Smoothing factor")
    parser.add_argument("-d", "--offdiag", dest="off_diag",action="store",required=False, type=int, default=2,
                        help="How many off diagonal needs to be included")
    parser.add_argument("-r", "--resfrag", dest="resfrag",action="store",required=True,
                        help="HiCpro resfrag file")
    parser.add_argument("-a", "--annotation", dest="valid_chroms",action="store",required=True,
                        help="Chromosomes annotation file")

    args = parser.parse_args()

    output_file = os.path.abspath(args.output_file)
    temporary_loc = os.path.dirname(output_file)



    CSR_mat,frag_index,frag_prop,frag_amount,valid_chroms,chroms_offsets = interaction_to_sparse.HiCpro_to_sparse(args.hicpro_results,args.resfrag,args.valid_chroms,temporary_loc,"bedgraph_temp")

    diagonal , num_reads = sparse_to_peaks.extract_diagonal(CSR_mat,args.off_diag)
    diagonal = sparse_to_peaks.moving_average(diagonal,args.smoothing)


    # combine that with the fragment prop and you will have it :)
    #checking length of fragprop and diagonal to be the same
    if len(frag_prop) != len(diagonal):
        raise Exception("something went wrong, diagonal size different")
    with open(output_file, "w") as bdg_file:
        for i in range(1,len(diagonal)-1):
            if frag_prop[i-1][0] != frag_prop[i+1][0]:
                continue
            bdg_file.write("{}\t{}\t{}\t{}\n".format(frag_prop[i-1][0],math.floor((frag_prop[i-1][2]+frag_prop[i-1][1])/2),math.floor((frag_prop[i][2]+frag_prop[i][1])/2),diagonal[i]))


if __name__=="__main__":
    main()