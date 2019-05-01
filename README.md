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



#### differential peak analysis

With the --keepdiff flag enabled it will produce a temporary file that can be used to do differential peak analysis
to integrate all samples together use the diff_peaks utility. this utility will look for all the correct files in a folder, merge the peaks at a fragment site level a produce a table with the signal in each peak from each sample. This can then be imported in R and analysed using DESeq2 or other differential expression analysis tools.







## Authors

## License

## Citation
