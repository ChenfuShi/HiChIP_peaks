#########################################
# Author: Chenfu Shi
# Email: chenfu.shi@postgrad.manchester.ac.uk


# collects data from single runs and compiles in a table easily analyzable with DESeq2
#########################################

def main():
    import os
    import argparse
    import logging
    import datetime
    import pickle
    import numpy
    try:
        #this works only when installed
        from hichip_tool import interaction_to_sparse,sparse_to_peaks
    except:
        import interaction_to_sparse 
        import sparse_to_peaks

    parser = argparse.ArgumentParser(description="input directory with outputfiles from peak_call and create table for differential analysis. Make sure to activate --keep_diff in the previous step!")

    parser.add_argument("-i", "--input", dest="hichip_tool_results",action="store",required=True,
                        help="directory containing previous step results")
    parser.add_argument("-o", "--output", dest="output_directory",action="store",required=True,
                        help="Output file, will write temporary files in that directory")

    parser.add_argument("-r", "--resfrag", dest="resfrag",action="store",required=True,
                        help="HiCpro resfrag file")          
    parser.add_argument("-a", "--annotation", dest="sizes",action="store",required=False, default=None,
                        help="HiCpro chromosome annotation file, default uses human chromosomes, excludes chrY")

    parser.add_argument("-m", "--minimum", dest="minimum",action="store",required=True,type=int,default=1,
                        help="How many samples need to be peak to be considered peak for analysis")          



    args = parser.parse_args()


    hichip_tool_results = os.path.abspath(args.hichip_tool_results)
    resfrag = os.path.abspath(args.resfrag)
    sizes = args.sizes
    output_dir = os.path.abspath(args.output_directory)

    minimum_coverage = args.minimum

    if not os.path.isdir(hichip_tool_results):
        raise Exception("couldn't find folder containing previous results")
    if not os.path.isfile(resfrag):
        raise Exception("annotation files couldn't be opened")
    if sizes!=None and not os.path.isfile(sizes) :
        raise Exception("annotation files couldn't be opened")

    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    # prepare variables

    frag_index,frag_prop,frag_amount,valid_chroms, chroms_offsets = interaction_to_sparse.Read_resfrag(resfrag,sizes)

    # look for files in the input directory. this will specifically look for files ending with diffpeak_data.pickle
    input_files = [fn for fn in os.listdir(hichip_tool_results) if fn.endswith("diffpeak_data.pickle")]
    if len(input_files) < 2:
        raise Exception("Couldn't find at least 2 files for differential expression analysis in the input directory!")
    names = [i.replace("diffpeak_data.pickle","") for i in input_files]
  
    input_smoothed_diagonals = []
    input_backgrounds = []
    input_refined_peaks = []

    for current_file in input_files:
        with open(current_file, "rb") as current_pickle:
            current_smoothed_diagonal, current_refined_peaks, current_background=pickle.load(current_pickle)
        input_smoothed_diagonals.append(current_smoothed_diagonal.copy())
        input_backgrounds.append(current_background.copy())
        input_refined_peaks.append(current_refined_peaks.copy())

    # merge the peaks
    stacked_data = numpy.stack(input_refined_peaks)
    stacked_peaks=stacked_data.sum(axis=0)
    merged_peaks=numpy.where(stacked_peaks >= minimum_coverage, 1, 0)

    # clean the peaks
    # clean up peak calling by removing peaks that are only 1 width
    peaks_string = "".join(["1" if x else "0" for x in merged_peaks])
    cleaned_string = peaks_string.replace("00100","00000")
    cleaned_string = cleaned_string.replace("11011","11111")
    # return list
    merged_peaks=[int(x) for x in list(cleaned_string)]


if __name__=="__main__":
    main()