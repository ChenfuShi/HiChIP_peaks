# HiChIP peaks

This package can be used to find enriched peak regions from HiChIP datasets that can then be used as an input to available loop calling tools or to do differential peak analysis.

It takes the HiC-Pro output and converts it to a restriction site level resolution map. It then selects reads within a specified number of restriction sites from the diagonal(default = 2) and models the background as a negative binomial. It calls peaks regions that significantly exceed the background.
The output is a list of peaks with their properties and a bedgraph at a restriction site level resolution that describes the reads per site.
Using the differential analysis command it can be used to create a consensus peakset and then identify differentially bound regions across samples.

Results from this package can then be used for further analysis and as a peaks dataset input for various loop calling software.

### Table of Contents 
- [HiChIP peaks](#HiChIP-peaks)
    - [Table of Contents](#Table-of-Contents)
  - [Getting started](#Getting-started)
    - [Installation](#Installation)
    - [Usage](#Usage)
      - [Peak calling](#Peak-calling)
      - [Example run](#Example-run)
      - [Differential peak analysis](#Differential-peak-analysis)
  - [Authors](#Authors)
  - [License](#License)
  - [Citation](#Citation)

## Getting started

### Installation

The package requires bedtools to run. The package can then be installed through pip

```
pip install hichip-peaks
```

We suggest using conda environments to avoid cluttering

```
conda create --name hichip-peaks python=3.7 bedtools pip
pip install hichip-peaks
```


### Usage
#### Peak calling
Clean the raw reads and align using HiC-Pro with normal settings making sure that these settings are set as follows(for MboI digested libraries):

```
#######################################################################
## Digestion Hi-C
#######################################################################

GENOME_FRAGMENT = MboI_resfrag_hg38.bed
LIGATION_SITE = GATCGATC
MIN_FRAG_SIZE = 
MAX_FRAG_SIZE =
MIN_INSERT_SIZE =
MAX_INSERT_SIZE =

#######################################################################
## Hi-C processing
#######################################################################

MIN_CIS_DIST =
GET_ALL_INTERACTION_CLASSES = 1
GET_PROCESS_SAM = 0
RM_SINGLETON = 1
RM_MULTI = 1
RM_DUP = 1
```

Use the peak_call command

```
usage: peak_call [-h] -i HICPRO_RESULTS -o OUTPUT_DIRECTORY -r RESFRAG
                 [-p PREFIX] [-f FDR] [-a SIZES] [-t TEMPORARY_LOC]
                 [-w THREADS] [-k] [-d] [-s OFF_DIAG] [-x] [-c]

Peak calling from HiChIP data

optional arguments:
  -h, --help            show this help message and exit
  -i HICPRO_RESULTS, --input HICPRO_RESULTS
                        HiC-Pro results directory containing validPairs file
                        and others
  -o OUTPUT_DIRECTORY, --output OUTPUT_DIRECTORY
                        Output directory
  -r RESFRAG, --resfrag RESFRAG
                        HiCpro resfrag file
  -p PREFIX, --prefix PREFIX
                        Output file name prefix, if not provided will be name
                        of HiC-Pro results directory
  -f FDR, --FDR FDR     False discovery rate, default = 0.01
  -a SIZES, --annotation SIZES
                        HiCpro chromosome annotation file, default uses human
                        chromosomes, excludes chrY
  -t TEMPORARY_LOC, --temporary_loc TEMPORARY_LOC
                        Temporary directory. If not supplied will be output
                        directory
  -w THREADS, --worker_threads THREADS
                        Number of threads, minimum 4. Warning: Increasing this
                        significantly increases RAM usage
  -k, --keep_temp       Keep temporary files
  -d, --keep_diff       Prepare files for differential analysis
  -s OFF_DIAG, --offdiag OFF_DIAG
                        How many off diagonal needs to be included (default =
                        2)
  -x, --chromX          Want to compensate Sex chromosomes weights? Requires
                        specify annotation(SIZES) containing chrX and chrY
  -c, --class_store     Store sparse site_matrix object for further use

```

This command needs the HICPRO_RESULTS/hic_results/data/sample/ output folder where all the valid pairs files are.
The command requires that all the files in that folder are present, including the .REPairs, SCPairs and DEPairs file.

**IMPORTANT** Please make sure that there is no other file named *allValidPairs* in the folder outside of the one containing the pairs, as this will cause the software to fail with *list index out of range* error. This includes a file that previous versions of HiC-pro created in the directory (allValidPairs.mergestat)

This command will generate the following files:

- log.log file, contains all the inputs used, the logs and the quality metrics calculated such as number of peaks called and fraction of reads in peaks.
- bedgraph.bdg file, containing the coverage track with all the reads used in the peak calling step.
- peaks.bed file, contains all the peaks called. The 3 extra information columns are:
  * average signal in the peak
  * maximum signal in the peak
  * -log10 of p value of the peak
- report.pdf, contains some useful plots and quality metrics.

If enabled this command will also generate:
- diffpeak_data.pickle, file containing the information necessary for differential peak analysis.
- mat_obj.pickle, file containing the sparse matrix representation at a restriction site level of all the interactions. Currently in development but you can look into site_matrix_class.py to understand how it works.

#### Example run
Assuming the data is in HICPRO_RESULTS/hic_results/data/sample/ and installation of HiC-Pro is in HICPRO_dir/
You can run the software with the following command:

```
peak_call -i HICPRO_RESULTS/hic_results/data/sample/ -o ./results -r HICPRO_dir/annotation/MboI_resfrag_hg38.bed 

```
Some example results can be found in [example_results](https://github.com/ChenfuShi/HiChIP_peaks/tree/master/example_results).

#### Differential peak analysis

Run the previous commands with the --keepdiff flag enabled. This will produce a temporary file that can be used with the diff_peaks command to integrate all samples together. This utility will look for all the correct files in a specified folder, merge the peaks at a fragment site level a produce a table with the signal in each peak from each sample. This can then be imported in R or others and analysed using DESeq2 or other differential expression analysis tools. See example R script for inspiration.

```
usage: diff_peaks [-h] -i hichip_peaks_RESULTS -o OUTPUT_FILE -r RESFRAG
                  [-a SIZES] [-m MINIMUM]

input directory with outputfiles from peak_call and create table for
differential analysis. Make sure to activate --keep_diff in the previous step!

optional arguments:
  -h, --help            show this help message and exit
  -i hichip_peaks_RESULTS, --input hichip_peaks_RESULTS
                        directory containing previous step results
  -o OUTPUT_FILE, --output OUTPUT_FILE
                        Output file
  -r RESFRAG, --resfrag RESFRAG
                        HiCpro resfrag file
  -a SIZES, --annotation SIZES
                        HiCpro chromosome annotation file, default uses human
                        chromosomes, excludes chrY
  -m MINIMUM, --minimum MINIMUM
                        How many samples need to be peak to be considered peak
                        for analysis
```


## Authors

This package was developed by Chenfu Shi<sup>1</sup>, Magnus Rattray<sup>2,3</sup> and Gisela Orozco<sup>1,3</sup> at the University of Manchester.

1) Centre for Genetics and Genomics Versus Arthritis. Division of Musculoskeletal and Dermatolog-ical Sciences, School of Biological Sciences, Faculty of Biology, Medicine and Health, The Universi-ty of Manchester, UK.
2) Division of Informatics, Imaging and Data Sciences, Faculty of Biology, Medicine and Health, University of Manchester, UK.
3) NIHR Manchester Biomedical Research Centre, Manchester University NHS Foundation Trust, Manchester Academic Health Science Centre, Manchester, UK.


This work was funded by the Wellcome Trust (award references 207491/Z/17/Z and 215207/Z/19/Z), Versus Arthritis (award reference 21754), NIHR Manchester BRC and the Medical Research Council (award reference MR/N00017X/1).



## License

The software is released with a BSD-3-Clause License

```
BSD-3-Clause License
Copyright 2019 Chenfu Shi
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
```

## Citation

Please cite our paper when using this package for your research!
```
Shi,C. et al. (2019) HiChIP-Peaks: A HiChIP peak calling algorithm. bioRxiv, 682781.
https://doi.org/10.1101/682781
```
