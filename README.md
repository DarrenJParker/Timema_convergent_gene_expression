# Timema_convergent_gene_expression

This is the repository for the collected scripts used in the study *"Repeated evolution of asexuality involves convergent gene expression changes."* currently under review.

## Components

### DATA

* Contains read counts, GO terms, dN/dS and pN/pS estimates for input to the scripts below. 

### Differential expression analyses

* **10sp_edgeR.R** | Script to identify gene expression changes between sexual and asexual species, using 10 species orthologs as a reference
* **cross_sp_edgeR.R** | Script to identify gene expression changes between sexual and asexual species, using each species as a reference

### GO term analyses

* **10sp_topGO.R** | Script to perform GO-term enrichment analyses using 10 species orthologs
* **cross_sp_topGO.R** | Script to perform GO-term enrichment analyses using each species as a reference

### pN/pS and dN/dS analyses

* **10sp_pNpSdNdS.R** | Script to analyse pN/pS and dN/dS for convergent vs background genes

### Number of convergent genes permutation analysis

* **readcount_sex_asex_randomiser.py** | Script to randomly switches the assignment of reproductive mode (sexual or asexual) within a species-pair for the readcount file (./Data/readcounts/10sp_orth_readcounts.csv)
* 

### Additional scripts

* **B2G_to_topGO.py** | Script for converting Blast2GO output into a format usable by topGO
* **Get_GO_term_parent_and_child_overlap_adjuster.py** | Script for dealing with topographically close GO terms 
* **Get_GO_term_parent_and_child_overlap_adjuster_test_data_expl.R** | Script explaining the logic of Get_GO_term_parent_and_child_overlap_adjuster.py
* **super_exact_test_multitest_corrector.py** | Script to correct the p-values of SuperExactTest for multiple tests
* **super_exact_test_table_parser.py** | Script to tidy up the output of super_exact_test_multitest_corrector.py


## Infomation on running scripts

* All scripts should be run from the directory they are in. Output directories will be created to store output as the code is run. 
* All python scripts were made using python 3.5. All contain help information which can be displayed by specifying no command line arguments.
* R code for cross_sp_edgeR.R and cross_sp_topGO.R are expected to be run in a loop for each reference species. For example (using bash):

```
for i in Tbi Tte Tce Tms Tcm Tsi Tpa Tge Tps Tdi; do
    Rscript cross_sp_edgeR.R $i
done
```

### for running the Number of convergent genes permutation analysis

* First make the datasets

```
python readcount_sex_asex_randomiser.py -i Data/readcounts/10sp_orth_readcounts.csv -N 10000 -o rand
```

* Run Expression analysis on each file 

```
for file in ./rand_ASRAND/*.csv; do
    Rscript 10sp_EdgeR_for_randomised_datasets.R $file rand_out
done
```

* collect up the number of convergent genes

```
python Nconvergentgenes_out_tidier.py -i rand_out_Nconvergentgenes_out/ -o rand
```


* Plot histograms to see distribution using 10sp_EdgeR_for_randomised_datasets_plots.R





## Abbreviations

### Species names:

Species name | Abbreviation | Reproductive mode 
--- | --- | --- 
*Timema bartmani* | Tbi | sexual 
*Timema tahoe* | Tte | asexual
*Timema cristinae* | Tce | sexual 
*Timema monikensis* | Tms | asexual
*Timema poppensis* | Tps | sexual 
*Timema douglasi* | Tdi | asexual
*Timema californicum* | Tcm | sexual 
*Timema shepardi* | Tsi | asexual
*Timema podura* | Tpa | sexual 
*Timema genevievae* | Tge | asexual

### Tissues:

Tissue | Abbreviation 
--- | --- 
Whole-body| WB
Reproductive tract | RT
Legs | LG




