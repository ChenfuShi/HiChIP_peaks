# domain_caller_site

This tool can be used to find enriched peak regions from HiChIP datasets.
It takes the HiC-Pro output and converts it to a restriction site level resolution map. Then it models the enrichment as a negative binomial and finds regions that exceed the background levels.

It will output a list of peaks with their properties and a bedgraph at a restriction site level resolution that describes the reads per site.

with the --keepdiff flag enabled it will produce a temporary file that can be used to do differential peak analysis
to integrate all samples together use the diff_peaks utility. this utility will look for all the correct files in a folder, merge the peaks at a fragment site level a produce a table with the signal in each peak from each sample. This can then be imported in R and analysed using DESeq2 or other differential expression analysis tools.