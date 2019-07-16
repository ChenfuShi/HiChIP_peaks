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
    import pandas
    import math
    try:
        #this works only when installed
        from hichip_peaks import interaction_to_sparse,sparse_to_peaks
    except:
        import interaction_to_sparse 
        import sparse_to_peaks

    parser = argparse.ArgumentParser(description="input directory with outputfiles from peak_call and create table for differential analysis. Make sure to activate --keep_diff in the previous step!")

    parser.add_argument("-i", "--input", dest="hichip_peaks_results",action="store",required=True,
                        help="directory containing previous step results")
    parser.add_argument("-o", "--output", dest="output_file",action="store",required=True,
                        help="Output file")
    parser.add_argument("-r", "--resfrag", dest="resfrag",action="store",required=True,
                        help="HiCpro resfrag file")          
    parser.add_argument("-a", "--annotation", dest="sizes",action="store",required=False, default=None,
                        help="HiCpro chromosome annotation file, default uses human chromosomes, excludes chrY")
    parser.add_argument("-m", "--minimum", dest="minimum",action="store",required=False,type=int,default=1,
                        help="How many samples need to be peak to be considered peak for analysis")          

    args = parser.parse_args()

    hichip_peaks_results = os.path.abspath(args.hichip_peaks_results)
    resfrag = os.path.abspath(args.resfrag)
    sizes = args.sizes
    output_file = os.path.abspath(args.output_file)

    minimum_coverage = args.minimum

    if not os.path.isdir(hichip_peaks_results):
        raise Exception("couldn't find folder containing previous results")
    if not os.path.isfile(resfrag):
        raise Exception("annotation files couldn't be opened")
    if sizes!=None and not os.path.isfile(sizes) :
        raise Exception("annotation files couldn't be opened")

    # prepare variables

    frag_index,frag_prop,frag_amount,valid_chroms, chroms_offsets = interaction_to_sparse.Read_resfrag(resfrag,sizes)

    # look for files in the input directory. this will specifically look for files ending with diffpeak_data.pickle
    input_files = [fn for fn in os.listdir(hichip_peaks_results) if fn.endswith("diffpeak_data.pickle")]
    if len(input_files) < 2:
        raise Exception("Couldn't find at least 2 files for differential expression analysis in the input directory!")
    names = [i.replace("diffpeak_data.pickle","") for i in input_files]
    print("Found {} files which will be merged".format(len(names)))
    input_smoothed_diagonals = []
    input_backgrounds = []
    input_refined_peaks = []

    # load one file at a time and append info into the lists
    for current_file in input_files:
        with open(os.path.join(hichip_peaks_results,current_file), "rb") as current_pickle:
            current_smoothed_diagonal, current_refined_peaks, current_background=pickle.load(current_pickle)
        input_smoothed_diagonals.append(current_smoothed_diagonal.copy())
        input_backgrounds.append(current_background.copy())
        input_refined_peaks.append(current_refined_peaks.copy())

    # stack the data
    stacked_diagonals = numpy.stack(input_smoothed_diagonals)
    stacked_backgrounds = numpy.stack(input_backgrounds)
    stacked_signal = numpy.subtract(stacked_diagonals,stacked_backgrounds)

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

    #get tuple of peaks. filtering at least minimum size of 3 restriction sites
    peaks=[]

    i = 0
    while i < len(merged_peaks):
        if merged_peaks[i] == 1:
            start = i
            while merged_peaks[i] == 1:
                i=i+1
            end = i
            if end - start > 2 :
                peaks.append((start,end))
        i=i+1

    print("Merged {} peaks from samples".format(len(peaks)))

    peaks_dict = {}
    for i in range(len(peaks)):
        peak_id = "peak_" + str(i)
        peak_location = frag_prop[peaks[i][0]][0] + ":" + str(math.floor((frag_prop[peaks[i][0]-1][1]+frag_prop[peaks[i][0]-1][2])/2))+"-"+str(math.floor((frag_prop[peaks[i][1]][1]+frag_prop[peaks[i][1]][2])/2))
        weigths = []
        for sample in range(len(names)):
            weigths.append(stacked_signal[sample][peaks[i][0]:peaks[i][1]].sum())
        peaks_dict[peak_id] = [peak_location] + weigths

    peaks_df=pandas.DataFrame.from_dict(peaks_dict, orient="index")
    peaks_df.columns = ["peak_location"] + names

    print("writing the peaks weights to {}".format(output_file))
    peaks_df.to_csv(output_file)


if __name__=="__main__":
    main()