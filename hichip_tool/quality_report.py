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

def quality_report(np_p_vals, smoothed_diagonal, output_dir, prefix):
    """takes all useful data and generates a good report"""

    pdf = matplotlib.backends.backend_pdf.PdfPages(os.path.join(output_dir, prefix + "report.pdf"))





    pdf.close()
    return None