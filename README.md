# HoldFst
 **H**eterzyg**o**sity scores for **l**ow-coverage **D**NA with Fixation Index (**Fst**)
 \
 \
HoldFst is tool for analyzing low-coverage whole genome sequencing. Without genotype information, standard Genome Wide Association Studies (GWAS) cannot be performed on low-coverage sequencing. HoldFst provides an alternative method for determining regions associated with a phenotype. Pooled low-coverage sequencing can immulate high-coverage sequencing so that heterozygosity and fixation can be calculated in windows across the genome. Regions of high significance can then be investigated for biological relevence. After narrowing down regions of significance using **HoldFst**, specific genes can be investigated in the region that could have biological implications related to the phenotype investigated.

## Quick Installation of HoldFst
First, install devtools (for installing GitHub packages) if it isn't already installed:
```
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
```
Then, install HoldFst. 
To include the vignette in the package installation add `build_vignettes = T`. (Note: this will increase install time)
```
devtools::install_github("aanderson54/HoldFst/HoldFst")
```
Lastly, load HoldFst
```
library(HoldFst)
```




## Usage
HoldFst uses a single VCF as the input file. VCFs can be downloaded from publicly available data or generated from new data. The VCF should have standard VCF columns as shown below. The columns named *population1* and *population2* should be replaced with the names of the two populations that you wish to compare as they appear in the vcf. 
```
CHROM POS ID  REF ALT QUAL  FILTER  INFO  FORMAT  *Population1* *Population2*
```
The package comes with a small vcf sample that contains variant calls from Chromsome1 from a population of striped and spotted Akagera zebras. To load the sample data,
```
data(Chromosome1)
```
To run HoldFst,
```
results<-HoldFst(file,P,Q,bin,plot)
```
* **file** is the file path to the vcf (ex. file="~/sample.vcf")
* **P** is the name of the normal population
* **Q** is the name of the variant population
* **bin** is the size of the bins that should be calculated across the genome (default=500,000)
* **plot** logical. Should the results be plotted?

To get the most significant regions of the genome,
```
topwindows(results,"all")
```
* **results** the object that you save the HoldFst data to
* The options for **sig** are Fst, Phet, Qhet, Diff, and all



## Studies that have used Heterozygosity and Fixation Index
Freedman, A. H., Lohmueller, K. E., & Wayne, R. K. (2016). Evolutionary history, selective sweeps, and deleterious variation in the dog. Annual Review of Ecology, Evolution, and Systematics, 47, 73-96. [link](https://www.annualreviews.org/doi/abs/10.1146/annurev-ecolsys-121415-032155?casa_token=8cdV2_R4HmgAAAAA:-nF6I-DL0b_oOjOcpFCiV9ZyroKMubY_dKIXmo73J2YSS8yiyxKgKmqqiG97RzUs2K73Xhd4k6mc)

Robinson, J. A., Räikkönen, J., Vucetich, L. M., Vucetich, J. A., Peterson, R. O., Lohmueller, K. E., & Wayne, R. K. (2019). Genomic signatures of extensive inbreeding in Isle Royale wolves, a population on the threshold of extinction. Science advances, 5(5), eaau0757. [link](https://advances.sciencemag.org/content/5/5/eaau0757.abstract)

Zhi, D., Da, L., Liu, M., Cheng, C., Zhang, Y., Wang, X., ... & Long, X. (2018). Whole genome sequencing of Hulunbuir short-tailed sheep for identifying candidate genes related to the short-tail phenotype. G3: Genes, Genomes, Genetics, 8(2), 377-383.
[link](https://www.g3journal.org/content/8/2/377.abstract)
\
\
\
...and many more
