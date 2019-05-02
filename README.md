# domain_caller_site

This package can be used to find enriched peak regions from HiChIP datasets that can then be used as an input to available loop calling tools or to do differential peak analysis.

It takes the HiC-Pro output and converts it to a restriction site level resolution map. It then selects reads within a specified number of restriction sites from the diagonal(default = 2) and models the background as a negative binomial. It calls peaks regions that significantly exceed the background.
The output is a list of peaks with their properties and a bedgraph at a restriction site level resolution that describes the reads per site.
Using the differential analysis command it can be used to create a consensus peakset and then identify differentially bound regions across samples.


## Getting started

### Installation

This package can be installed through pip

```
pip install mypackage
```

Or through bioconda

```
conda install mypackage
```

We suggest using conda environments to avoid cluttering


### Usage
#### peak calling
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
                 [-w THREADS] [-k] [-d] [-s OFF_DIAG]

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
                        Number of threads, minimum 4
  -k, --keep_temp       Keep temporary files
  -d, --keep_diff       Prepare files for differential analysis
  -s OFF_DIAG, --offdiag OFF_DIAG
                        How many off diagonal needs to be included

```

This command needs the HiC-pro_results/hic_results/data/sample output folder where all the valid pairs files are.
The command requires that all the files in that folder are present, including the .REPairs, SCPairs and DEPairs file.




#### differential peak analysis

Run the previous commands with the --keepdiff flag enabled. This will produce a temporary file that can be used with the diff_peaks  command to integrate all samples together. This utility will look for all the correct files in a folder, merge the peaks at a fragment site level a produce a table with the signal in each peak from each sample. This can then be imported in R or others and analysed using DESeq2 or other differential expression analysis tools. See example R script for inspiration.

```
usage: diff_peaks [-h] -i HICHIP_TOOL_RESULTS -o OUTPUT_FILE -r RESFRAG
                  [-a SIZES] [-m MINIMUM]

input directory with outputfiles from peak_call and create table for
differential analysis. Make sure to activate --keep_diff in the previous step!

optional arguments:
  -h, --help            show this help message and exit
  -i HICHIP_TOOL_RESULTS, --input HICHIP_TOOL_RESULTS
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


## License

## Citation
