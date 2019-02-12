# Timema_convergent_gene_expression

This is the repository for the collected scripts used in the study:

>Parker, D. J., Bast, J., Jalvingh, K., Dumas, Z., Robinson-Rechavi, M., Schwander, T. 2019. [Repeated evolution of asexuality involves convergent gene expression changes](https://academic.oup.com/mbe/article/36/2/350/5184916?guestAccessKey=8e0f0601-f224-45c4-801b-92d1f8622bd2). *Molecular Biology and Evolution*. 36: 350–364

Scripts for the paper are archived at http://dx.doi.org/10.5281/zenodo.2025853.

# Components

## DATA

* Contains read counts, GO terms, dN/dS and pN/pS estimates for input to the scripts below. 

## SCRIPTS

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
* **10sp_EdgeR_for_randomised_datasets.R** | EdgeR script to identify convergent gene expression changes between sexual and asexual species for the randomised datasets
* **Nconvergentgenes_out_tidier.py** | Collects output from 10sp_EdgeR_for_randomised_datasets.R runs
* **10sp_EdgeR_for_randomised_datasets_plots.R** | plots histograms of the permutation analysis

### Ornstein-Uhlenbeck (OU) models

* **split_expression_file_for_OU.py** | This script takes the logCPM files from 10sp_edgeR.R and produces files that can be used to run Ornstein-Uhlenbeck (OU) models in the R package ouch.
* **OU_models.R** | Script to run the OU models
* **python OU_summarise_and_FDR.py** | Collects output from OU models corrects for multiple tests

### Additional scripts

* **B2G_to_topGO.py** | Script for converting Blast2GO output into a format usable by topGO
* **Get_GO_term_parent_and_child_overlap_adjuster.py** | Script for dealing with topographically close GO terms 
* **Get_GO_term_parent_and_child_overlap_adjuster_test_data_expl.R** | Script explaining the logic of Get_GO_term_parent_and_child_overlap_adjuster.py
* **super_exact_test_multitest_corrector.py** | Script to correct the p-values of SuperExactTest for multiple tests
* **super_exact_test_table_parser.py** | Script to tidy up the output of super_exact_test_multitest_corrector.py


# Infomation on running scripts

## General 

* All scripts should be run from the directory they are in. Output directories will be created to store output as the code is run. 
* All python scripts were made using python 3.5. All contain help information which can be displayed by specifying no command line arguments.
* R code for cross_sp_edgeR.R and cross_sp_topGO.R are expected to be run in a loop for each reference species. For example (using bash):

```
for i in Tbi Tte Tce Tms Tcm Tsi Tpa Tge Tps Tdi; do
    Rscript cross_sp_edgeR.R $i
done
```

## To run the Number of convergent genes permutation analysis

* First make the datasets

```
python3.5 readcount_sex_asex_randomiser.py -i Data/readcounts/10sp_orth_readcounts.csv -N 10000 -o rand
```

* Run Expression analysis on each file 

```
for file in ./rand_ASRAND/*.csv; do
    Rscript 10sp_EdgeR_for_randomised_datasets.R $file rand_out
done
```

* collect up the number of convergent genes

```
python3.5 Nconvergentgenes_out_tidier.py -i rand_out_Nconvergentgenes_out/ -o rand
```


* Plot histograms to see distribution using 10sp_EdgeR_for_randomised_datasets_plots.R


## To run the Ornstein-Uhlenbeck (OU) models

* First the differential expression analyses using 10sp_edgeR.R must be run to obtain mean expression datasets for each tissue. This data files will be in the (newly created) directory `Gene_exp_for_OU_models`

* Split the expression file into a form that can be used by ouch

```
python3.5 split_expression_file_for_OU.py
```

* run OU models on each gene

```
## run in a loop in bash
for tissue in WB RT LG ; do
    in_dir=`echo $tissue"_GE_for_OU"`
    out_dir=`echo $tissue"_OU_out"`
    
    for file in ./$in_dir/*.txt; do
        Rscript OU_models.R $file $out_dir
	done
done
```

* Collect output and correct for multiple tests

```
python3.5 OU_summarise_and_FDR.py -w WB_OU_out -r RT_OU_out -l LG_OU_out
```

# Abbreviations

## Species names:

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

## Tissues:

Tissue | Abbreviation 
--- | --- 
Whole-body| WB
Reproductive tract | RT
Legs | LG
