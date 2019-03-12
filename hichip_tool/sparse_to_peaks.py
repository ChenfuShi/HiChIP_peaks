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
import functools

def sparse_to_peaks(CSR_mat,frag_index,frag_prop,frag_amount,valid_chroms,chroms_offsets,output_dir,prefix,FDR=0.10,threads=4,keeptemp=False):
    """Wrapper function to call individual funcitons"""

    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)


    print("#######################################")
    print("Extracting pairs for ChIP peaks calling")

    diagonal , num_reads = extract_diagonal(CSR_mat,2)
    print("Number of reads used in peak calling: {}".format(num_reads))
    if num_reads < 20000000:
        print("WARNING: number of reads used for peak calling is very low. Consider doing more sequencing")
    smoothed_diagonal = numpy.rint(moving_integration(diagonal,3)).astype(int) #### changed to 3

    print("#######################################")
    print("Identifying high confidence peaks to remove them from background modelling")

    quick_peaks = quick_call(smoothed_diagonal)

    refined_peaks , peak_p_vals , peaks_q_vals= refined_call(smoothed_diagonal,quick_peaks,frag_prop,FDR,threads)

    print("#######################################")
    print("Writing peaks and bedgraph to output folder")

    output_bed = os.path.join(output_dir, prefix + "peaks.bed")
    output_bedgraph =  os.path.join(output_dir, prefix + "bedgraph.bdg")
    bed_printout(frag_prop,smoothed_diagonal,refined_peaks,peak_p_vals,output_bed,output_bedgraph)
    
    if keeptemp==True:
        with open(os.path.join(output_dir, prefix + "peaks_variables.pi"),"wb") as picklefile:
            pickle.dump([smoothed_diagonal, refined_peaks ,quick_peaks, peak_p_vals , peaks_q_vals],picklefile)    


    #peaks returned is just a list with 0 and 1. proper bed file is saved
    return smoothed_diagonal, refined_peaks ,quick_peaks, peak_p_vals , peaks_q_vals

def moving_integration (values, window):
    weights = numpy.repeat(1.0, window)
    sma = numpy.convolve(values, weights, 'same')
    return sma
def moving_average (values, window):
    weights = numpy.repeat(1.0, window)/window
    sma = numpy.convolve(values, weights, 'same')
    return sma

def get_range(frag_prop, index, distance):
    """finds the index ranges for a specified distance range, used to get the local average of noise"""
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
    diagonal = CSR_mat.diagonal()/2
    num_reads = sum(diagonal)
    if window == 0:
        return numpy.array(diagonal),num_reads
    for i in range(1,window+1):
        off_diagonal = CSR_mat.diagonal(k=i).tolist()
        num_reads += sum(off_diagonal)
        diagonal = [sum(x) for x in zip(diagonal, [0]*i + off_diagonal, off_diagonal + [0]*i)]
    return numpy.array(diagonal),num_reads

def quick_call(smoothed_diagonal):
    """calls the peaks using a very simple genomic average"""

    average_signal = numpy.mean(smoothed_diagonal)
    quick_p_vals=[]
    poisson_pre_pvals = [scipy.stats.poisson.sf(x, average_signal) + scipy.stats.poisson.pmf(x, average_signal) for x in range(numpy.max(smoothed_diagonal)+1)]

    for res_site in smoothed_diagonal.tolist(): 
        quick_p_vals.append(poisson_pre_pvals[res_site])

    #quick_peaks, correct_q_vals = statsmodels.stats.multitest.fdrcorrection(quick_p_vals, alpha = 0.01)

    quick_peaks = [True if x < 0.00000001 else False for x in quick_p_vals] #seems to work better
    
    # matplotlib.pyplot.plot([1000 if x else 0 for x in quick_peaks])
    # matplotlib.pyplot.plot(smoothed_diagonal)
    # matplotlib.pyplot.show()

    return quick_peaks

def parallel_expected_background(frag_prop,smoothed_diagonal,noise_filter,group_lengths,max_interpolated,min_interpolated,diagonal_mean,size_function,mean_mode_diff,i):
    """parallel version of expected background"""
    start_index , end_index = get_range(frag_prop,i,10000)
    local_background = get_local_background(noise_filter,smoothed_diagonal,start_index,end_index)
    local_length = group_lengths[i]
    if local_length > max_interpolated:
        local_length = max_interpolated-1
    if local_length < min_interpolated:
        local_length = min_interpolated+1
    if local_background > diagonal_mean:
        return size_function(local_length) + mean_mode_diff + local_background - diagonal_mean
    else:
        return size_function(local_length) + mean_mode_diff

def parallel_negative_binomial(expected_background,nb_n,smoothed_diagonal,site_index):
    nb_p = nb_n/(expected_background[site_index]+nb_n)
    return scipy.stats.nbinom.sf(smoothed_diagonal[site_index], nb_n,nb_p) + scipy.stats.nbinom.pmf(smoothed_diagonal[site_index] , nb_n,nb_p) 






def refined_call(smoothed_diagonal, quick_peaks, frag_prop,FDR,threads):
    """use previous peaks to refine model and then call peaks. creates a list with expected noise based on measures. poisson distribution won't work, need to increase variance.
    then clean up isolated stuff and return peaks"""

    print("#######################################")
    print("Model background noise as a negative binomial")

    lengths = [x[3] for x in frag_prop] 
    group_lengths = moving_integration(lengths, 2)  ###changed to 2, so the 2 fragments within
    min_allowed_size = math.floor(numpy.percentile(group_lengths, 1))
    max_allowed_size = math.floor(numpy.percentile(group_lengths, 99))
    # add 1s to quick peaks to exclude them from the noise modelling
    noise_filter = quick_peaks.copy()
    for i in range(len(group_lengths)):
        if group_lengths[i] > max_allowed_size or group_lengths[i] < min_allowed_size:
            noise_filter[i] = 1

    # select only the bits that give you noise. On which to model stuff
    noise_lengths = list(itertools.compress(group_lengths, [not i for i in noise_filter]))
    noise_diagonal = list(itertools.compress(smoothed_diagonal, [not i for i in noise_filter]))


    # estimate overdispersion parameter from data
    nbinom_data = statsmodels.api.NegativeBinomial(noise_diagonal,numpy.ones(len(noise_diagonal)),disp=False)
    nb = nbinom_data.fit()
    nb_const, nb_alpha = nb.params

    print("Negative binomial overdispersion parameter: {}".format(nb_alpha))
    print("#######################################")
    print("Identify effect of fragment size bias")

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
    #ynew = size_function(xnew)

    # matplotlib.pyplot.plot([size_function(x) for x in range(min_interpolated,max_interpolated)])
    # matplotlib.pyplot.show()

    print("#######################################")
    print("Estimating expected background levels from local background and fragment size")

    # associate a mean with every site, from the size distribution that is not a mean, it's a mode. maybe correct for the local mean as well? how do you actually associate both
    # calculate mean. if local mean is higher add that amount to the result of the interpolation
    # smaller than or bigger than just use the latest value

    size_mean = numpy.mean(predicted_diagonals)
    diagonal_mean = numpy.mean(noise_diagonal)

    mean_mode_diff = diagonal_mean - size_mean 

    # print("doing parallel stuff now")
    parallel_partial_expected_back = functools.partial(parallel_expected_background,frag_prop,smoothed_diagonal,noise_filter,group_lengths,max_interpolated,min_interpolated,diagonal_mean,size_function,mean_mode_diff)
    pool = multiprocessing.Pool(threads)
    expected_background = pool.map(parallel_partial_expected_back, range(len(smoothed_diagonal)),chunksize=50000)
    pool.close()
    pool.join()

    # print("parallel stuff finishes now")
    # expected_background=[]
    # for i in range(len(smoothed_diagonal)):
    #     start_index , end_index = get_range(frag_prop,i,10000)
    #     local_background = get_local_background(noise_filter,smoothed_diagonal,start_index,end_index)
    #     local_length = group_lengths[i]
    #     if local_length > max_interpolated:
    #         local_length = max_interpolated-1
    #     if local_length < min_interpolated:
    #         local_length = min_interpolated+1
    #     if local_background > diagonal_mean:
    #         expected_background.append(size_function(local_length) + mean_mode_diff + local_background - diagonal_mean)
    #     else:
    #         expected_background.append(size_function(local_length) + mean_mode_diff)
    # if expected_background==expected_background2:
    #     print("seems fine")


    # matplotlib.pyplot.plot(smoothed_diagonal)
    # matplotlib.pyplot.plot(expected_background)
    # matplotlib.pyplot.plot(expected_background2)
    # matplotlib.pyplot.show()

    print("#######################################")
    print("Identifying enriched regions using negative binomial model")

    # run peak calling using a negative binomial model, input the p and mean calculated using the mean and the dispersion parameter from the nb fit
    nb_p_vals = []
    nb_n = 1/nb_alpha
    # for site_index in range(len(smoothed_diagonal)): 
    #     nb_p = nb_n/(expected_background[site_index]+nb_n)
    #     nb_p_vals.append(scipy.stats.nbinom.sf(smoothed_diagonal[site_index], nb_n,nb_p) + scipy.stats.nbinom.pmf(smoothed_diagonal[site_index] , nb_n,nb_p) )

    # print("parallel negative binomial test")

    parallel_partial_nb_test = functools.partial(parallel_negative_binomial,expected_background,nb_n,smoothed_diagonal)
    pool = multiprocessing.Pool(threads)
    nb_p_vals = pool.map(parallel_partial_nb_test, range(len(smoothed_diagonal)),chunksize=50000)
    pool.close()
    pool.join()
    # if nb_p_vals2==nb_p_vals:
    #     print("seems fine")    


    matplotlib.pyplot.hist(nb_p_vals,bins=50)
    matplotlib.pyplot.show()

    # MAYBE false discovery rate depending on how bad it looks like or set a proper p value. macs uses different thing for FDR and its possible only if using controls.
    nb_peaks, nb_q_vals = statsmodels.stats.multitest.fdrcorrection(nb_p_vals, alpha = FDR)

    # nb_peaks = [True if x < 0.001 else False for x in nb_p_vals] #seems to work better

    # clean up peak calling by removing peaks that are only 1 width
    peaks_string = "".join(["1" if x else "0" for x in nb_peaks])
    cleaned_string = peaks_string.replace("00100","00000")
    cleaned_string = cleaned_string.replace("11011","11111")
    # return list
    refined_peaks=[int(x) for x in list(cleaned_string)]

    print("#######################################")
    print("Refined peak calling done")
    
    # matplotlib.pyplot.plot([1000 if x==1 else 0 for x in refined_peaks])
    # matplotlib.pyplot.plot(smoothed_diagonal)
    # matplotlib.pyplot.show()


    if len(expected_background) != len(smoothed_diagonal) or len(smoothed_diagonal) != len(group_lengths) or len(refined_peaks) != len(group_lengths):
        raise Exception("something happened in the lengths of the various vectors")
    return refined_peaks , nb_p_vals, nb_q_vals



def bed_printout(frag_prop,smoothed_diagonal,refined_peaks,peak_p_vals,output_bed,output_bedgraph):
    """print out a bed file with refined peaks, also add as a score the fold change of the highest point"""

    with open(output_bed + ".temp", "w") as output_file:
        for i in range(1,len(smoothed_diagonal)-1):
            if frag_prop[i-1][0] != frag_prop[i+1][0]:
                continue
            if refined_peaks[i] == 1:
                output_file.write("{}\t{}\t{}\t{}\t{:10.15f}\n".format(frag_prop[i-1][0],math.floor((frag_prop[i-1][2]+frag_prop[i-1][1])/2),math.floor((frag_prop[i][2]+frag_prop[i][1])/2),smoothed_diagonal[i],-math.log10(peak_p_vals[i])))
    bedmerge_command = "bedtools merge -i " + output_bed + ".temp -c 4,4,5 -o mean,max,max> " + output_bed 
    subprocess.check_call(bedmerge_command ,shell=True)
    os.remove(output_bed + ".temp")

    with open(output_bedgraph, "w") as bdg_file:
        for i in range(1,len(smoothed_diagonal)-1):
            if frag_prop[i-1][0] != frag_prop[i+1][0]:
                continue
            bdg_file.write("{}\t{}\t{}\t{}\n".format(frag_prop[i-1][0],math.floor((frag_prop[i-1][2]+frag_prop[i-1][1])/2),math.floor((frag_prop[i][2]+frag_prop[i][1])/2),smoothed_diagonal[i]))

                    



if __name__=="__main__":
    """test functions here"""
    import pickle
    CSR_mat = scipy.sparse.load_npz('../domain_caller_site/testdata/sparse_matrix_mumbach_non_reassigned_chr1.npz')
    with open("../domain_caller_site/testdata/variables.pi","rb") as picklefile:
        frag_index,frag_prop,frag_amount,valid_chroms,chroms_offsets = pickle.load(picklefile)
    output_dir = os.path.abspath("./testdata")
    smoothed_diagonal, refined_peaks ,quick_peaks, peak_p_vals , peaks_q_vals = sparse_to_peaks(CSR_mat,frag_index,frag_prop[:584662],frag_amount,valid_chroms,chroms_offsets,output_dir,"testdata",threads=6)

    with open("./testdata/peaks_chr1_mumbach.pi","wb") as picklefile:
        pickle.dump([smoothed_diagonal, refined_peaks ,quick_peaks, peak_p_vals , peaks_q_vals],picklefile)


