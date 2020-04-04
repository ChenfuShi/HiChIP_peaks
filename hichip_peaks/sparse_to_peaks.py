#########################################
# Author: Chenfu Shi
# Email: chenfu.shi@postgrad.manchester.ac.uk


# Function inputs a HiChIP CSR_mat and extracts the peaks. 
# Extracts short range interactions from these peaks, forms a vector. Calls significantly enriched fragments then cleans up and merges peaks together. 
# Returns list of peaks as fragment numbers

#########################################

import scipy
import scipy.sparse
import scipy.stats , scipy.interpolate
import statsmodels.stats.multitest , statsmodels , statsmodels.api
import math
import numpy
import os
import re
import multiprocessing
import subprocess
import matplotlib.pyplot
import itertools
#import functools
import logging
import pickle

def sparse_to_peaks(CSR_mat,frag_index,frag_prop,frag_amount,valid_chroms,chroms_offsets,output_dir,prefix,off_diag,chromX,FDR=0.01,threads=4,keeptemp=False):
    """Wrapper function to call individual funcitons
    smoothed_diagonal : (list) that contains the number of reads from the diagonal + n off_diags and then is smoothed using the moving integration
    refined_peaks : (list) peaks output from the negative binomial model
    quick_peaks : (list) peaks output from the poisson model
    peak_p_vals : (list) p values for each restriction site
    peaks_q_vals : (list) q, corrected p values for each restriction site
    expected_background: (list) expected background for each site calculated using negative binomial model and restriction fragment size bias
    """

    # create output directory
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    logging.info("#######################################")
    logging.info("Extracting pairs for ChIP peaks calling")

    # Extract diagonal from sparse matrix
    diagonal , num_reads = extract_diagonal(CSR_mat,off_diag) #off diag here
    
    logging.info("Number of reads used in peak calling: {}".format(num_reads/2))
    if num_reads < 30000000:
        logging.warning("WARNING: number of reads used for peak calling is very low. Consider doing more sequencing")
    
    # moving integration of diagonal
    smoothed_diagonal = numpy.rint(moving_integration(diagonal,((off_diag-1)*2)+1)).astype(int) #### changed to 3 smoothing factor, it is one less than the number of off sites

    # if check for sex chromosomes is active run check. Will return correct diagonal with doubled counts on X and Y if male sample recognized
    if chromX == True:
        smoothed_diagonal = checkX(smoothed_diagonal,chroms_offsets,frag_amount)

    logging.info("#######################################")
    logging.info("Identifying high confidence peaks to remove them from background modelling")

    # quick call of peaks using poisson model test with very stringent p value
    quick_peaks = quick_call(smoothed_diagonal)

    # slow call of peaks using negative binomial model
    refined_peaks , peak_p_vals , peaks_q_vals, expected_background= refined_call(smoothed_diagonal,quick_peaks,frag_prop,FDR,off_diag,threads)

    logging.info("#######################################")
    logging.info("Writing peaks and bedgraph to output folder")

    # Create output files
    output_bed = os.path.join(output_dir, prefix + "peaks.bed")
    output_bedgraph =  os.path.join(output_dir, prefix + "bedgraph.bdg")
    bed_printout(frag_prop,smoothed_diagonal,refined_peaks,peak_p_vals,output_bed,output_bedgraph,expected_background,keeptemp)
    
    if keeptemp==True:
        with open(os.path.join(output_dir, prefix + "peaks_variables.pi"),"wb") as picklefile:
            pickle.dump([smoothed_diagonal, refined_peaks ,quick_peaks, peak_p_vals , peaks_q_vals,expected_background],picklefile)    

    # peaks returned is just a list with 0 and 1. proper bed file is saved
    return smoothed_diagonal, refined_peaks ,quick_peaks, peak_p_vals , peaks_q_vals ,expected_background


def moving_integration (values, window):
    weights = numpy.repeat(1.0, window)
    sma = numpy.convolve(values, weights, 'same')
    return sma
def moving_average (values, window):
    weights = numpy.repeat(1.0, window)/window
    sma = numpy.convolve(values, weights, 'same')
    return sma


def checkX(smoothed_diagonal,chroms_offsets,frag_amount):
    """This function checks for the presence of the chromosome X and Y and checks if Y has more than 0.2 counts. 
    if it has counts then it doubles counts from X and Y chroms"""
    chroms_prop = {}
    for i in chroms_offsets.keys(): 
        chroms_prop[i] = (chroms_offsets[i], chroms_offsets[i] + frag_amount[i])
    if "chrX" in chroms_prop.keys() and "chrY" in chroms_prop.keys():
        if smoothed_diagonal[chroms_prop["chrY"][0]:chroms_prop["chrY"][1]].mean() > smoothed_diagonal.mean() * 0.2:
            logging.info("Recognized male sample. Doubling counts from X and Y chromosomes")
            new_diagonal = []   
            for i in chroms_prop.keys(): 
                if i == "chrY" or i == "chrX":
                    new_diagonal.extend(list(map(lambda x :x*2, smoothed_diagonal[chroms_prop[i][0]:chroms_prop[i][1]])))
                else:
                    new_diagonal.extend(smoothed_diagonal[chroms_prop[i][0]:chroms_prop[i][1]])
            return new_diagonal
        else:
            logging.info("Recognized female sample. Do nothing")
            return smoothed_diagonal

    elif "X" in chroms_prop.keys() and "Y" in chroms_prop.keys():
        if smoothed_diagonal[chroms_prop["Y"][0]:chroms_prop["Y"][1]].mean() > smoothed_diagonal.mean() * 0.2:
            logging.info("Recognized male sample. Doubling counts from X and Y chromosomes")  
            new_diagonal = []   
            for i in chroms_prop.keys(): 
                if i == "Y" or i == "X":
                    new_diagonal.extend(list(map(lambda x :x*2, smoothed_diagonal[chroms_prop[i][0]:chroms_prop[i][1]])))
                else:
                    new_diagonal.extend(smoothed_diagonal[chroms_prop[i][0]:chroms_prop[i][1]])
            return new_diagonal
        else:
            logging.info("Recognized female sample. Do nothing")
            return smoothed_diagonal
    
    else:
        logging.info("No chromosome X or Y in annotation file. Do nothing")
        return smoothed_diagonal


def get_range(frag_prop, index, distance):
    """THIS FUNCTION IS NOT USED finds the index ranges for a specified distance range, used to get the local average of noise"""
    chromosome = frag_prop[index][0]
    start = index
    end = index
    initial_pos = frag_prop[index][1]
    max_end = len(frag_prop)
    while True:
        start = start - 1
        if start < 0:
            start = 0
            break
        if initial_pos - frag_prop[start][1] > distance or chromosome != frag_prop[start][0]:
            start = start + 1
            break
    while True:
        end = end + 1
        if end >= max_end:
            end = max_end-1
            break
        if frag_prop[end][1] - initial_pos > distance or chromosome != frag_prop[end][0]:
            end = end - 1
            break

    return start, end


def get_local_background(signal_list, smoothed_diagonal, start_index, end_index):
    """THIS FUNCTION IS NOT USED gets local background around the start and end index"""
    background = 0
    used_sites = 0
    for i in range(start_index,end_index+1):
        if signal_list[i] == 0:
            background += smoothed_diagonal[i]
            used_sites += 1
    if used_sites == 0:
        return 0
    local_background = background/used_sites
    return local_background


def extract_diagonal(CSR_mat,window):
    """extract the diagonal including the sum of the window in all directions"""
    diagonal = CSR_mat.diagonal()
    num_reads = sum(diagonal)
    if window == 0:
        return numpy.array(diagonal),num_reads
    for i in range(1,window+1):
        off_diagonal = CSR_mat.diagonal(k=i).tolist()
        num_reads += sum(off_diagonal)*2
        diagonal = [sum(x) for x in zip(diagonal, [0]*i + off_diagonal, off_diagonal + [0]*i)]
    return numpy.array(diagonal),num_reads


def quick_call(smoothed_diagonal):
    """calls the peaks using a very simple genomic average and a poisson test"""

    # calculate genomic average
    average_signal = numpy.mean(smoothed_diagonal)
    
    # Pre-calculate all the p values for every possible signal intensity. This speeds up the calculation massively
    poisson_pre_pvals = [scipy.stats.poisson.sf(x, average_signal) + scipy.stats.poisson.pmf(x, average_signal) for x in range(numpy.max(smoothed_diagonal)+1)]
    
    # Assign p values for all sites
    quick_p_vals=[]
    for res_site in smoothed_diagonal.tolist(): 
        quick_p_vals.append(poisson_pre_pvals[res_site])

    # Filter sites by p value 10^-8
    quick_peaks = [True if x < 0.00000001 else False for x in quick_p_vals]
    
    return quick_peaks


# parallel implementation of expected background estimation. Initialize thread with data and then worker calls the function
expback_aux_data = None
def initializer_parallel_expected_background(group_lengths,max_interpolated,min_interpolated,diagonal_mean,size_function,mean_mode_diff):
    global expback_aux_data
    expback_aux_data = [group_lengths,max_interpolated,min_interpolated,diagonal_mean,size_function,mean_mode_diff]
def worker_parallel_expected_background(i):
    global expback_aux_data
    assert expback_aux_data is not None
    return parallel_expected_background(expback_aux_data[0],expback_aux_data[1],expback_aux_data[2],expback_aux_data[3],expback_aux_data[4],expback_aux_data[5],i)
def parallel_expected_background(group_lengths,max_interpolated,min_interpolated,diagonal_mean,size_function,mean_mode_diff,i):
    """parallel version of expected background. Uses the size_function to estimate how many reads do you expect because of the size of the fragments"""
    #start_index , end_index = get_range(frag_prop,i,100000)
    #local_background = get_local_background(noise_filter,smoothed_diagonal,start_index,end_index) ##seems that whatever i do this always make it worse
    local_background = 0
    local_length = group_lengths[i]
    if local_length > max_interpolated:
        local_length = max_interpolated-1
    if local_length < min_interpolated:
        local_length = min_interpolated+1
    if local_background > diagonal_mean: # this is a fake if with the current implementation
        return size_function(local_length) + mean_mode_diff #+ local_background - diagonal_mean
    else:
        return size_function(local_length) + mean_mode_diff


# parallel implementation of negative binomial test. Initialize thread with data and then worker calls the function
nb_aux_data = None
def initializer_parallel_negative_binomial(expected_background,nb_n,smoothed_diagonal):
    global nb_aux_data
    nb_aux_data = [expected_background,nb_n,smoothed_diagonal]
def worker_parallel_negative_binomial(site_index):
    global nb_aux_data
    assert nb_aux_data is not None
    return parallel_negative_binomial(nb_aux_data[0],nb_aux_data[1],nb_aux_data[2],site_index)
def parallel_negative_binomial(expected_background,nb_n,smoothed_diagonal,site_index):
    nb_p = nb_n/(expected_background[site_index]+nb_n)
    return scipy.stats.nbinom.sf(smoothed_diagonal[site_index], nb_n,nb_p) + scipy.stats.nbinom.pmf(smoothed_diagonal[site_index] , nb_n,nb_p) 


def refined_call(smoothed_diagonal, quick_peaks, frag_prop,FDR,off_diag,threads):
    """use previous peaks to refine model and then call peaks. creates a list with expected noise based on measures. poisson distribution won't work, need to increase variance.
    then clean up isolated stuff and return peaks"""

    logging.info("#######################################")
    logging.info("Model background noise as a negative binomial")

    # get the lengths of the fragments and set the size of the group of fragments within the moving integration of the smoothing
    lengths = [x[3] for x in frag_prop] 
    group_lengths = moving_integration(lengths, (off_diag-1)*2)  ###changed to 2, so the 2 fragments within the smoothing factor, calculate it from the window used to do the smoothing
    # setting boundaries for lowess estimation
    min_allowed_size = math.floor(numpy.percentile(group_lengths, 1))
    max_allowed_size = math.floor(numpy.percentile(group_lengths, 99))
    
    #  set 1s to quick peaks to exclude them from the background modelling
    noise_filter = quick_peaks.copy()
    for i in range(len(group_lengths)):
        if group_lengths[i] > max_allowed_size or group_lengths[i] < min_allowed_size:
            noise_filter[i] = 1

    # select only the bits that give you background on which to model stuff
    noise_lengths = list(itertools.compress(group_lengths, [not i for i in noise_filter]))
    noise_diagonal = list(itertools.compress(smoothed_diagonal, [not i for i in noise_filter]))

    # estimate overdispersion parameter from data. nb_const is basically mean
    nbinom_data = statsmodels.api.NegativeBinomial(noise_diagonal,numpy.ones(len(noise_diagonal)),disp=False)
    nb = nbinom_data.fit()
    nb_const, nb_alpha = nb.params

    logging.info("Negative binomial overdispersion parameter: {}".format(nb_alpha))
    logging.info("#######################################")
    logging.info("Identify effect of fragment size bias")

    # lowess fit the size distribution
    # subset of 200k fragments
    idx = numpy.random.choice(len(noise_lengths), size=200000, replace=False)
    subset_lengths = [noise_lengths[n] for n in idx]
    subset_diagonal = [noise_diagonal[n] for n in idx]
    # find lowess prediction
    predicted_stuff = statsmodels.api.nonparametric.lowess(subset_diagonal, subset_lengths,return_sorted=True, frac=0.4 , delta=3.0 ) 
    predicted_lengths = list(zip(*predicted_stuff))[0]
    predicted_diagonals = list(zip(*predicted_stuff))[1]
    # use that prediction to interpolate a function to predict new data
    size_function = scipy.interpolate.interp1d(predicted_lengths, predicted_diagonals, bounds_error=True) 
    max_interpolated = max(predicted_lengths)
    min_interpolated = min(predicted_lengths)

    logging.info("#######################################")
    logging.info("Estimating expected background levels from fragment size")

    # from the size distribution that is not a mean, it's a mode. need to calculate the difference so that it corrects when using the size function. 

    size_mean = numpy.mean(predicted_diagonals)
    diagonal_mean = numpy.mean(noise_diagonal)

    mean_mode_diff = diagonal_mean - size_mean 
    
    # parallel assignment of expected background. initializer sends the lists and then map uses just an index
    pool = multiprocessing.Pool(threads,initializer_parallel_expected_background,[group_lengths,max_interpolated,min_interpolated,diagonal_mean,size_function,mean_mode_diff])
    expected_background = pool.map(worker_parallel_expected_background, range(len(smoothed_diagonal)))
    pool.close()
    pool.join()

    logging.info("#######################################")
    logging.info("Identifying enriched regions using negative binomial model")

    # run peak calling using a negative binomial model, input the p and mean calculated using the mean and the dispersion parameter from the nb fit
    nb_n = 1/nb_alpha
    # parallel calculation of p value from negative binomial. initializer function sends lists and then map uses just an index
    pool = multiprocessing.Pool(threads,initializer_parallel_negative_binomial,[expected_background,nb_n,smoothed_diagonal])
    nb_p_vals = pool.map(worker_parallel_negative_binomial, range(len(smoothed_diagonal)))
    pool.close()
    pool.join()

    # false discovery rate correction
    nb_peaks, nb_q_vals = statsmodels.stats.multitest.fdrcorrection(nb_p_vals, alpha = FDR)

    # clean up peak calling by removing peaks that are only 1 width and filling holes that are only 1 width
    peaks_string = "".join(["1" if x else "0" for x in nb_peaks])
    cleaned_string = peaks_string.replace("00100","00000")
    cleaned_string = cleaned_string.replace("11011","11111")
    # return list
    refined_peaks=[int(x) for x in list(cleaned_string)]

    logging.info("#######################################")
    logging.info("Refined peak calling done")
    
    # quick check that everything is alright
    if len(expected_background) != len(smoothed_diagonal) or len(smoothed_diagonal) != len(group_lengths) or len(refined_peaks) != len(group_lengths):
        raise Exception("something happened in the lengths of the various vectors")
    return refined_peaks , nb_p_vals, nb_q_vals, expected_background


def bed_printout(frag_prop,smoothed_diagonal,refined_peaks,peak_p_vals,output_bed,output_bedgraph,expected_background,keeptemp):
    """print out a bed file with refined peaks, also add as a score the fold change of the highest point
    also prints out a bedgraph with the signal used for the peak calling"""
    # prints all fragments and then uses bedtools to merge the fragments that are contigous.
    with open(output_bed + ".temp", "w") as output_file:
        for i in range(1,len(smoothed_diagonal)-1):
            if frag_prop[i-1][0] != frag_prop[i+1][0]:
                continue
            if refined_peaks[i] == 1:
                output_file.write("{}\t{}\t{}\t{}\t{:10.15f}\n".format(frag_prop[i-1][0],math.floor((frag_prop[i-1][2]+frag_prop[i-1][1])/2),math.floor((frag_prop[i][2]+frag_prop[i][1])/2),max(smoothed_diagonal[i]-expected_background[i],0),-math.log10(peak_p_vals[i])))
    bedmerge_command = "bedtools merge -i " + output_bed + ".temp -c 4,4,5 -o mean,max,max > " + output_bed 
    subprocess.check_call(bedmerge_command ,shell=True)
    if keeptemp == False:
        os.remove(output_bed + ".temp")

    # standard bedgraph format
    with open(output_bedgraph, "w") as bdg_file:
        for i in range(1,len(smoothed_diagonal)-1):
            if frag_prop[i-1][0] != frag_prop[i+1][0]:
                continue
            bdg_file.write("{}\t{}\t{}\t{}\n".format(frag_prop[i-1][0],math.floor((frag_prop[i-1][2]+frag_prop[i-1][1])/2),math.floor((frag_prop[i][2]+frag_prop[i][1])/2),smoothed_diagonal[i]))

                    



if __name__=="__main__":
    """test functions here"""
    import pickle
    logging.basicConfig(
        level=logging.INFO,
        format="%(levelname)s - %(message)s",
        handlers=[
        logging.StreamHandler()
    ]
    )
    CSR_mat = scipy.sparse.load_npz('./testdata/sparse_matrix.npz')
    with open("./testdata/variables.pi","rb") as picklefile:
        frag_index,frag_prop,frag_amount,valid_chroms,chroms_offsets = pickle.load(picklefile)
    output_dir = os.path.abspath("./testdata")
    smoothed_diagonal, refined_peaks ,quick_peaks, peak_p_vals , peaks_q_vals ,expected_background= sparse_to_peaks(CSR_mat,frag_index,frag_prop,frag_amount,valid_chroms,chroms_offsets,output_dir,"testdata",2,threads=6)


    with open("./testdata/peaks_chr1_mumbach.pi","wb") as picklefile:
        pickle.dump([smoothed_diagonal, refined_peaks ,quick_peaks, peak_p_vals , peaks_q_vals],picklefile)


