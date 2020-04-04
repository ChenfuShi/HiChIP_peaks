#########################################
# Author: Chenfu Shi
# Email: chenfu.shi@postgrad.manchester.ac.uk


# prepare quality report using data
# number of reads in peaks, how can you calculate this?
# ip efficiency doesn't need a separate file
# p value distribution

# from log. take number of reads imported and number of reads used in peak calling
# number of peaks
#########################################

import os, sys
import matplotlib.pyplot
import matplotlib.backends.backend_pdf
import scipy
import scipy.sparse, scipy.stats
import numpy
import logging

def quality_report(peak_p_vals,refined_peaks, smoothed_diagonal, output_dir, prefix):
    """takes all useful data and generates a good report"""
    #compatibility with no X display
    matplotlib.pyplot.switch_backend('agg')
    
    pdf = matplotlib.backends.backend_pdf.PdfPages(os.path.join(output_dir, prefix + "report.pdf"))
    #basic numbers info
    fig, ax = matplotlib.pyplot.subplots(figsize=(9,5))
    ax.text(0.5, 0.95, 'Thank you for using my software',
            horizontalalignment='center',
            verticalalignment='center',
        fontsize=18)
    num_reads = round(sum(smoothed_diagonal)/6)
    ax.text(0.5, 0.7, 'Number of reads used for peak calling : {}'.format(num_reads),
            horizontalalignment='center',
            verticalalignment='center',
        fontsize=16)

    num_lines = sum(1 for line in open(os.path.join(output_dir, prefix + "peaks.bed")))
    ax.text(0.5, 0.5, 'Number of peaks called : {}'.format(num_lines),
            horizontalalignment='center',
            verticalalignment='center',
        fontsize=16)

    num_reads_in_peaks=round(sum([x for x,y in zip(smoothed_diagonal, refined_peaks) if y==1])/6)

    ax.text(0.5, 0.3, 'Number of reads in peaks : {} ({:.2%})'.format(num_reads_in_peaks,num_reads_in_peaks/num_reads),
            horizontalalignment='center',
            verticalalignment='center',
        fontsize=16)

    ax.axis("off")
    pdf.savefig(fig)

    logging.info('Number of reads used for peak calling : {}'.format(num_reads))
    logging.info('Number of peaks called : {}'.format(num_lines))
    logging.info('Number of reads in peaks : {} ({:.2%})'.format(num_reads_in_peaks,num_reads_in_peaks/num_reads))

    #p value distribution
    fig2, ax2 = matplotlib.pyplot.subplots(figsize=(9,5))
    ax2.spines["top"].set_visible(False)  
    ax2.spines["right"].set_visible(False)  
    matplotlib.pyplot.title("Negative binomial test p-value distribution \nCheck for uniform distribution", fontsize=17)    
    matplotlib.pyplot.ylabel("frequency", fontsize=14)
    matplotlib.pyplot.xlabel("p-values", fontsize=14)
    ax2.tick_params(labelsize=13)
    ax2.hist(peak_p_vals, bins=20, range=(0,1),density=True,color="#3F5D7D")

    pdf.savefig(fig2)


    pdf.close()
    return None