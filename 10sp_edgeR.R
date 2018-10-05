#### 10sp_edgeR.R

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
### libraries

library("edgeR")
library("VennDiagram")
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(cowplot)
library(stringr)
library(gtable)
library(pheatmap)
library(RColorBrewer)
require(vegan)
library(pvclust)
library(raster)
library("SuperExactTest")

print (sessionInfo())

# # R version 3.4.1 (2017-06-30)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: macOS High Sierra 10.13.3

# Matrix products: default
# BLAS: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRblas.0.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRlapack.dylib

# locale:
# [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

# attached base packages:
# [1] grid      stats     graphics  grDevices utils     datasets  methods   base     

# other attached packages:
 # [1] SuperExactTest_0.99.4 raster_2.5-8          sp_1.2-5              pvclust_2.0-0        
 # [5] vegan_2.4-4           permute_0.9-4         RColorBrewer_1.1-2    stringr_1.2.0        
 # [9] VennDiagram_1.6.17    futile.logger_1.4.3   edgeR_3.18.1          limma_3.32.10        
# [13] gtable_0.2.0          lattice_0.20-35       gridExtra_2.3         cowplot_0.8.0        
# [17] pheatmap_1.0.8        ggplot2_2.2.1        

# loaded via a namespace (and not attached):
 # [1] Rcpp_0.12.13         cluster_2.0.6        magrittr_1.5         MASS_7.3-47         
 # [5] munsell_0.4.3        colorspace_1.3-2     rlang_0.1.6          plyr_1.8.4          
 # [9] tools_3.4.1          parallel_3.4.1       nlme_3.1-131         mgcv_1.8-22         
# [13] lambda.r_1.2         lazyeval_0.2.1       tibble_1.3.4         Matrix_1.2-11       
# [17] futile.options_1.0.0 stringi_1.1.5        compiler_3.4.1       scales_0.5.0        
# [21] locfit_1.5-9.1  

Intersect <- function (x) {  
  # Multiple set version of intersect
  # x is a list
  if (length(x) == 1) {
    unlist(x)
  } else if (length(x) == 2) {
    intersect(x[[1]], x[[2]])
  } else if (length(x) > 2){
    intersect(x[[1]], Intersect(x[-1]))
  }
}


#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
# Data

rawdata_10sp <- read.csv("DATA/readcounts/10sp_orth_readcounts.csv", check.names=FALSE, stringsAsFactors=FALSE)
head(rawdata_10sp)
colnames(rawdata_10sp)

##### output

dir.create("10sp_DE_output")

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
## into DGE structure - all species together first
## Tissues seperatly 
## JUST Females
## Note tissue-specific - ALL samples are mated


#### Whole-body samples
y_WB_10sp_UF <- DGEList(counts=rawdata_10sp[,c(
"Tbi_SF_WB_Re1","Tbi_SF_WB_Re2","Tbi_SF_WB_Re3","Tte_AF_WB_Re1","Tte_AF_WB_Re2","Tte_AF_WB_Re3",
"Tce_SF_WB_Re1","Tce_SF_WB_Re2","Tce_SF_WB_Re3","Tms_AF_WB_Re1","Tms_AF_WB_Re2","Tms_AF_WB_Re3",
"Tcm_SF_WB_Re1","Tcm_SF_WB_Re2","Tcm_SF_WB_Re3","Tsi_AF_WB_Re1","Tsi_AF_WB_Re2","Tsi_AF_WB_Re3",
"Tpa_SF_WB_Re1","Tpa_SF_WB_Re2","Tpa_SF_WB_Re3","Tge_AF_WB_Re1","Tge_AF_WB_Re2","Tge_AF_WB_Re3",
"Tps_SF_WB_Re1","Tps_SF_WB_Re2","Tps_SF_WB_Re3","Tdi_AF_WB_Re1","Tdi_AF_WB_Re2","Tdi_AF_WB_Re3"
)], genes=rawdata_10sp[,1:1])


#### Reproductive tract samples
y_RT_10sp_UF <- DGEList(counts=rawdata_10sp[,c(
"Tbi_SF_RT_Re1","Tbi_SF_RT_Re2","Tbi_SF_RT_Re3","Tte_AF_RT_Re1","Tte_AF_RT_Re2","Tte_AF_RT_Re3",
"Tce_SF_RT_Re1","Tce_SF_RT_Re2","Tce_SF_RT_Re3","Tms_AF_RT_Re1","Tms_AF_RT_Re2","Tms_AF_RT_Re3",
"Tcm_SF_RT_Re1","Tcm_SF_RT_Re2","Tcm_SF_RT_Re3","Tsi_AF_RT_Re1","Tsi_AF_RT_Re2","Tsi_AF_RT_Re3",
"Tpa_SF_RT_Re1","Tpa_SF_RT_Re2","Tpa_SF_RT_Re3","Tge_AF_RT_Re1","Tge_AF_RT_Re2","Tge_AF_RT_Re3",
"Tps_SF_RT_Re1","Tps_SF_RT_Re2","Tps_SF_RT_Re3","Tdi_AF_RT_Re1","Tdi_AF_RT_Re2","Tdi_AF_RT_Re3"
)], genes=rawdata_10sp[,1:1])

#### Leg samples
y_LG_10sp_UF <- DGEList(counts=rawdata_10sp[,c(
"Tbi_SF_LG_Re1","Tbi_SF_LG_Re2","Tbi_SF_LG_Re3","Tte_AF_LG_Re1","Tte_AF_LG_Re2","Tte_AF_LG_Re3",
"Tce_SF_LG_Re1","Tce_SF_LG_Re2","Tce_SF_LG_Re3","Tms_AF_LG_Re1","Tms_AF_LG_Re2","Tms_AF_LG_Re3",
"Tcm_SF_LG_Re1","Tcm_SF_LG_Re2","Tcm_SF_LG_Re3","Tsi_AF_LG_Re1","Tsi_AF_LG_Re2","Tsi_AF_LG_Re3",
"Tpa_SF_LG_Re1","Tpa_SF_LG_Re2","Tpa_SF_LG_Re3","Tge_AF_LG_Re1","Tge_AF_LG_Re2","Tge_AF_LG_Re3",
"Tps_SF_LG_Re1","Tps_SF_LG_Re2","Tps_SF_LG_Re3","Tdi_AF_LG_Re1","Tdi_AF_LG_Re2","Tdi_AF_LG_Re3"
)], genes=rawdata_10sp[,1:1])



#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### filter and normalisation (code)

filt_and_norm_maj <- function(y,cpm_cut,cut_in_Nsams){
	
	cat("\nNumber number of genes / samples in orig data\n")
	print(dim(y)) ### number of genes / samples
	head(cpm(y)) 
	keep <- 
	  rowSums(cpm(y[,1:3])> cpm_cut) >= cut_in_Nsams & 
      rowSums(cpm(y[,4:6])> cpm_cut) >= cut_in_Nsams &
      rowSums(cpm(y[,7:9])> cpm_cut) >= cut_in_Nsams &    
      rowSums(cpm(y[,10:12])> cpm_cut) >= cut_in_Nsams &     
      rowSums(cpm(y[,13:15])> cpm_cut) >= cut_in_Nsams &    
      rowSums(cpm(y[,16:18])> cpm_cut) >= cut_in_Nsams &      
      rowSums(cpm(y[,19:21])> cpm_cut) >= cut_in_Nsams &    
      rowSums(cpm(y[,22:24])> cpm_cut) >= cut_in_Nsams &     
      rowSums(cpm(y[,25:27])> cpm_cut) >= cut_in_Nsams &
      rowSums(cpm(y[,28:30])> cpm_cut) >= cut_in_Nsams
   

	y <- y[keep,]

	y$samples$lib.size <- colSums(y$counts) 
	y <- calcNormFactors(y)
	cat("\nLib norm factors\n")
	print(y$samples)	
	cat("\nNumber number of genes / samples after filtering\n")
	print(dim(y))
	
	return(y)

}


#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### filter and normalisation (samples)

y_WB_10sp <- filt_and_norm_maj(y_WB_10sp_UF,0.5,2)
# Number number of genes / samples after filtering
# 2985   30

y_RT_10sp <- filt_and_norm_maj(y_RT_10sp_UF,0.5,2)
# Number number of genes / samples after filtering
# 2754   30

y_LG_10sp <- filt_and_norm_maj(y_LG_10sp_UF,0.5,2)
# Number number of genes / samples after filtering
#  2737   30


#################################################################################################
### Calculate mean logCPMs for each sample for OU analysis

CPM_WB_10sp <- as.data.frame(cpm(y_WB_10sp, log = T))
rownames(CPM_WB_10sp) <- as.character(y_WB_10sp$genes$genes)

CPM_RT_10sp <- as.data.frame(cpm(y_RT_10sp, log = T))
rownames(CPM_RT_10sp) <- as.character(y_RT_10sp$genes$genes)

CPM_LG_10sp <- as.data.frame(cpm(y_LG_10sp, log = T))
rownames(CPM_LG_10sp) <- as.character(y_LG_10sp$genes$genes)


mean_CPM_calc <- function(df,vec_of_col_names,base_out){
	
	### make matrix
	
	a1 <- vec_of_col_names[1]
	
	print(a1)
	
	Rep1 <- eval(parse(text=paste(df,"$",vec_of_col_names[1],sep="")))
	Rep2 <- eval(parse(text=paste(df,'$',vec_of_col_names[2],sep='')))
	Rep3 <- eval(parse(text=paste(df,'$',vec_of_col_names[3],sep='')))

	mat_0 <- as.matrix(cbind(Rep1, Rep2, Rep3))
	df <- eval(parse(text=df))
	
	mean_list <- c() 
	
	for(i in 1:nrow(df)) {   
		a <- mat_0[i,]
		d <- mean(a) 
		mean_list <- c(mean_list,d)
	}
	

	df4 <- as.data.frame(cbind(mean_list))

	colnames(df4) <- paste(base_out, "CPM_mean", sep = "_")	
	rownames(df4) <- rownames(df)

	return(df4)
}


### WB 
WB_10sp_mean_CPM <- cbind(
mean_CPM_calc("CPM_WB_10sp", c("Tbi_SF_WB_Re1", "Tbi_SF_WB_Re2", "Tbi_SF_WB_Re3"), "Tbi_SF_WB"),
mean_CPM_calc("CPM_WB_10sp", c("Tte_AF_WB_Re1", "Tte_AF_WB_Re2", "Tte_AF_WB_Re3"), "Tte_AF_WB"),
mean_CPM_calc("CPM_WB_10sp", c("Tce_SF_WB_Re1", "Tce_SF_WB_Re2", "Tce_SF_WB_Re3"), "Tce_SF_WB"),
mean_CPM_calc("CPM_WB_10sp", c("Tms_AF_WB_Re1", "Tms_AF_WB_Re2", "Tms_AF_WB_Re3"), "Tms_AF_WB"),
mean_CPM_calc("CPM_WB_10sp", c("Tcm_SF_WB_Re1", "Tcm_SF_WB_Re2", "Tcm_SF_WB_Re3"), "Tcm_SF_WB"),
mean_CPM_calc("CPM_WB_10sp", c("Tsi_AF_WB_Re1", "Tsi_AF_WB_Re2", "Tsi_AF_WB_Re3"), "Tsi_AF_WB"),
mean_CPM_calc("CPM_WB_10sp", c("Tpa_SF_WB_Re1", "Tpa_SF_WB_Re2", "Tpa_SF_WB_Re3"), "Tpa_SF_WB"),
mean_CPM_calc("CPM_WB_10sp", c("Tge_AF_WB_Re1", "Tge_AF_WB_Re2", "Tge_AF_WB_Re3"), "Tge_AF_WB"),
mean_CPM_calc("CPM_WB_10sp", c("Tps_SF_WB_Re1", "Tps_SF_WB_Re2", "Tps_SF_WB_Re3"), "Tps_SF_WB"),
mean_CPM_calc("CPM_WB_10sp", c("Tdi_AF_WB_Re1", "Tdi_AF_WB_Re2", "Tdi_AF_WB_Re3"), "Tdi_AF_WB"))


### RT 
RT_10sp_mean_CPM <- cbind(
mean_CPM_calc("CPM_RT_10sp", c("Tbi_SF_RT_Re1", "Tbi_SF_RT_Re2", "Tbi_SF_RT_Re3"), "Tbi_SF_RT"),
mean_CPM_calc("CPM_RT_10sp", c("Tte_AF_RT_Re1", "Tte_AF_RT_Re2", "Tte_AF_RT_Re3"), "Tte_AF_RT"),
mean_CPM_calc("CPM_RT_10sp", c("Tce_SF_RT_Re1", "Tce_SF_RT_Re2", "Tce_SF_RT_Re3"), "Tce_SF_RT"),
mean_CPM_calc("CPM_RT_10sp", c("Tms_AF_RT_Re1", "Tms_AF_RT_Re2", "Tms_AF_RT_Re3"), "Tms_AF_RT"),
mean_CPM_calc("CPM_RT_10sp", c("Tcm_SF_RT_Re1", "Tcm_SF_RT_Re2", "Tcm_SF_RT_Re3"), "Tcm_SF_RT"),
mean_CPM_calc("CPM_RT_10sp", c("Tsi_AF_RT_Re1", "Tsi_AF_RT_Re2", "Tsi_AF_RT_Re3"), "Tsi_AF_RT"),
mean_CPM_calc("CPM_RT_10sp", c("Tpa_SF_RT_Re1", "Tpa_SF_RT_Re2", "Tpa_SF_RT_Re3"), "Tpa_SF_RT"),
mean_CPM_calc("CPM_RT_10sp", c("Tge_AF_RT_Re1", "Tge_AF_RT_Re2", "Tge_AF_RT_Re3"), "Tge_AF_RT"),
mean_CPM_calc("CPM_RT_10sp", c("Tps_SF_RT_Re1", "Tps_SF_RT_Re2", "Tps_SF_RT_Re3"), "Tps_SF_RT"),
mean_CPM_calc("CPM_RT_10sp", c("Tdi_AF_RT_Re1", "Tdi_AF_RT_Re2", "Tdi_AF_RT_Re3"), "Tdi_AF_RT"))


### LG 
LG_10sp_mean_CPM <- cbind(
mean_CPM_calc("CPM_LG_10sp", c("Tbi_SF_LG_Re1", "Tbi_SF_LG_Re2", "Tbi_SF_LG_Re3"), "Tbi_SF_LG"),
mean_CPM_calc("CPM_LG_10sp", c("Tte_AF_LG_Re1", "Tte_AF_LG_Re2", "Tte_AF_LG_Re3"), "Tte_AF_LG"),
mean_CPM_calc("CPM_LG_10sp", c("Tce_SF_LG_Re1", "Tce_SF_LG_Re2", "Tce_SF_LG_Re3"), "Tce_SF_LG"),
mean_CPM_calc("CPM_LG_10sp", c("Tms_AF_LG_Re1", "Tms_AF_LG_Re2", "Tms_AF_LG_Re3"), "Tms_AF_LG"),
mean_CPM_calc("CPM_LG_10sp", c("Tcm_SF_LG_Re1", "Tcm_SF_LG_Re2", "Tcm_SF_LG_Re3"), "Tcm_SF_LG"),
mean_CPM_calc("CPM_LG_10sp", c("Tsi_AF_LG_Re1", "Tsi_AF_LG_Re2", "Tsi_AF_LG_Re3"), "Tsi_AF_LG"),
mean_CPM_calc("CPM_LG_10sp", c("Tpa_SF_LG_Re1", "Tpa_SF_LG_Re2", "Tpa_SF_LG_Re3"), "Tpa_SF_LG"),
mean_CPM_calc("CPM_LG_10sp", c("Tge_AF_LG_Re1", "Tge_AF_LG_Re2", "Tge_AF_LG_Re3"), "Tge_AF_LG"),
mean_CPM_calc("CPM_LG_10sp", c("Tps_SF_LG_Re1", "Tps_SF_LG_Re2", "Tps_SF_LG_Re3"), "Tps_SF_LG"),
mean_CPM_calc("CPM_LG_10sp", c("Tdi_AF_LG_Re1", "Tdi_AF_LG_Re2", "Tdi_AF_LG_Re3"), "Tdi_AF_LG"))


### export 

dir.create("Gene_exp_for_OU_models")

write.csv(WB_10sp_mean_CPM, file="Gene_exp_for_OU_models/WB_10sp_mean_CPM.csv", row.names=TRUE)
write.csv(RT_10sp_mean_CPM, file="Gene_exp_for_OU_models/RT_10sp_mean_CPM.csv", row.names=TRUE)
write.csv(LG_10sp_mean_CPM, file="Gene_exp_for_OU_models/LG_10sp_mean_CPM.csv", row.names=TRUE)


#################################################################################################
### get DE_gene_table (code)

## returns full table and sig DE genes as a vector
get_DE_genes <- function(fita,FDRa){
	TT1 = topTags(fita, n =3000000000)
	TT2 = TT1$table
	temp_sig <- subset(TT2, TT2$FDR <= FDRa)	
	sig_genes = temp_sig$genes
	sig_logFC = temp_sig$logFC
	
	N_sig_genes <- length(sig_genes)
	cat("Number of sig genes: ", N_sig_genes )
	
	r_list <- list("table" = TT2, "S_gene_list" = sig_genes, "S_logFC_list" = sig_logFC )
	return(r_list)
}


####### design matrix
sp_tiss = factor(c(
"Tbi","Tbi","Tbi","Tbi","Tbi","Tbi",
"Tce","Tce","Tce","Tce","Tce","Tce",
"Tcm","Tcm","Tcm","Tcm","Tcm","Tcm",
"Tpa","Tpa","Tpa","Tpa","Tpa","Tpa",
"Tps","Tps","Tps","Tps","Tps","Tps"
))

prop_tiss = factor(c(
"S","S","S","A","A","A",
"S","S","S","A","A","A",
"S","S","S","A","A","A",
"S","S","S","A","A","A",
"S","S","S","A","A","A"
))

data.frame(Sample=colnames(y_RT_10sp),sp_tiss,prop_tiss)

### design matrix WB
design_WB_10sp <- model.matrix(~sp_tiss * prop_tiss )
rownames(design_WB_10sp) <- colnames(y_WB_10sp)
design_WB_10sp

### design matrix RT
design_RT_10sp <- model.matrix(~sp_tiss * prop_tiss )
rownames(design_RT_10sp) <- colnames(y_RT_10sp)
design_RT_10sp

### design matrix LG
design_LG_10sp <- model.matrix(~sp_tiss * prop_tiss )
rownames(design_LG_10sp) <- colnames(y_LG_10sp)
design_LG_10sp


###### Est dispersion (samples)

y_WB_10sp <- estimateDisp(y_WB_10sp, design_WB_10sp )
y_RT_10sp <- estimateDisp(y_RT_10sp, design_RT_10sp)
y_LG_10sp <- estimateDisp(y_LG_10sp, design_LG_10sp)


plotBCV(y_WB_10sp)
y_WB_10sp$common.dispersion
sqrt(y_WB_10sp$common.dispersion) ## this is the (common) BCV

plotBCV(y_RT_10sp )
y_RT_10sp$common.dispersion
sqrt(y_RT_10sp$common.dispersion)

plotBCV(y_LG_10sp )
y_LG_10sp$common.dispersion
sqrt(y_LG_10sp$common.dispersion)


#####################################################################################################################################################################
###### fit glm 

fit_4inter_WB_10sp <- glmFit(y_WB_10sp, design_WB_10sp,robust=TRUE)
fit_4inter_RT_10sp <- glmFit(y_RT_10sp, design_RT_10sp,robust=TRUE)
fit_4inter_LG_10sp <- glmFit(y_LG_10sp, design_LG_10sp,robust=TRUE)

#### overall interaction

fit_4inter_WB_10sp_c <- glmLRT(fit_4inter_WB_10sp,coef=7:10)
fit_4inter_RT_10sp_c <- glmLRT(fit_4inter_RT_10sp,coef=7:10)
fit_4inter_LG_10sp_c <- glmLRT(fit_4inter_LG_10sp,coef=7:10)

TTT_WB_10sp_S_A_INTER <- get_DE_genes(fit_4inter_WB_10sp_c, 0.05)
TTT_RT_10sp_S_A_INTER <- get_DE_genes(fit_4inter_RT_10sp_c, 0.05)
TTT_LG_10sp_S_A_INTER <- get_DE_genes(fit_4inter_LG_10sp_c, 0.05)

####### export full tables

write.table(TTT_WB_10sp_S_A_INTER$table, "10sp_DE_output/TTT_WB_10sp_S_A_INTER.csv", sep = ',', quote = FALSE, row.names = FALSE)
write.table(TTT_RT_10sp_S_A_INTER$table,            "10sp_DE_output/TTT_RT_10sp_S_A_INTER.csv", sep = ',', quote = FALSE, row.names = FALSE)
write.table(TTT_LG_10sp_S_A_INTER$table,            "10sp_DE_output/TTT_LG_10sp_S_A_INTER.csv", sep = ',', quote = FALSE, row.names = FALSE)

##### get overall effects using additive model

### design matrix WB_match_tiss
design_WB_10sp_a <- model.matrix(~sp_tiss + prop_tiss )
rownames(design_WB_10sp_a) <- colnames(y_WB_10sp)
design_WB_10sp_a

### design matrix RT
design_RT_10sp_a <- model.matrix(~sp_tiss + prop_tiss )
rownames(design_RT_10sp_a) <- colnames(y_RT_10sp)
design_RT_10sp_a

### design matrix LG
design_LG_10sp_a <- model.matrix(~sp_tiss + prop_tiss )
rownames(design_LG_10sp_a) <- colnames(y_LG_10sp)
design_LG_10sp_a

fit_4oveeff_WB_10sp <- glmFit(y_WB_10sp, design_WB_10sp_a,robust=TRUE)
fit_4oveeff_RT_10sp <- glmFit(y_RT_10sp, design_RT_10sp_a,robust=TRUE)
fit_4oveeff_LG_10sp <- glmFit(y_LG_10sp, design_LG_10sp_a,robust=TRUE)

fit_4propeff_WB_10sp <- glmLRT(fit_4oveeff_WB_10sp,coef=6)
fit_4propeff_RT_10sp <- glmLRT(fit_4oveeff_RT_10sp,coef=6)
fit_4propeff_LG_10sp <- glmLRT(fit_4oveeff_LG_10sp,coef=6)

TTT_WB_10sp_S_A_propeff <- get_DE_genes(fit_4propeff_WB_10sp, 0.05)
TTT_RT_10sp_S_A_propeff <- get_DE_genes(fit_4propeff_RT_10sp, 0.05)
TTT_LG_10sp_S_A_propeff <- get_DE_genes(fit_4propeff_LG_10sp, 0.05)


####### export full tables

write.table(TTT_WB_10sp_S_A_propeff$table, "10sp_DE_output/TTT_WB_10sp_S_A_propeff.csv", sep = ',', quote = FALSE, row.names = FALSE)
write.table(TTT_RT_10sp_S_A_propeff$table, "10sp_DE_output/TTT_RT_10sp_S_A_propeff.csv", sep = ',', quote = FALSE, row.names = FALSE)
write.table(TTT_LG_10sp_S_A_propeff$table, "10sp_DE_output/TTT_LG_10sp_S_A_propeff.csv", sep = ',', quote = FALSE, row.names = FALSE)


#####################################################################################################################################################################
###### fit glm for contrasts (code) 

make_Group_matrix <- function(y, fact1,fact2){
	a_samp = data.frame(Sample=colnames(y),fact1,fact2)
	Group <- factor(paste(a_samp$fact1,a_samp$fact2, sep="."))
	print(Group)
	cbind(a_samp,Group=Group)
	cat(Group)
	G_design <- model.matrix(~0+Group)
	colnames(G_design) <- levels(Group)
	print(G_design)
	return(G_design)
}

glmy_fit <- function(y,G_design){
	fit_samples <- glmFit(y, G_design,robust=TRUE)
	return(fit_samples)
}


####### fit glm (samples)

mat_WB_10sp <- make_Group_matrix(y_WB_10sp,sp_tiss,prop_tiss)
fit_WB_10sp <- glmy_fit(y_WB_10sp,mat_WB_10sp)

mat_RT_10sp <- make_Group_matrix(y_RT_10sp,sp_tiss,prop_tiss)
fit_RT_10sp <- glmy_fit(y_RT_10sp,mat_RT_10sp)

mat_LG_10sp <- make_Group_matrix(y_LG_10sp,sp_tiss,prop_tiss)
fit_LG_10sp <- glmy_fit(y_LG_10sp,mat_LG_10sp)

#################################################################################################
### make contrasts


WB_10sp.contrasts <- makeContrasts(
Tbi_S_A	=	Tbi.A - Tbi.S ,   ### +ve FC = higher in ASEX 
Tce_S_A	=	Tce.A - Tce.S ,
Tcm_S_A	=	Tcm.A - Tcm.S ,
Tpa_S_A	=	Tpa.A - Tpa.S ,
Tps_S_A	=	Tps.A - Tps.S ,
levels=mat_WB_10sp)

Tbi_lrt_WB_10sp_S_A  <- glmLRT(fit_WB_10sp, contrast=WB_10sp.contrasts[,"Tbi_S_A"])
Tce_lrt_WB_10sp_S_A  <- glmLRT(fit_WB_10sp, contrast=WB_10sp.contrasts[,"Tce_S_A"])
Tcm_lrt_WB_10sp_S_A  <- glmLRT(fit_WB_10sp, contrast=WB_10sp.contrasts[,"Tcm_S_A"])
Tpa_lrt_WB_10sp_S_A  <- glmLRT(fit_WB_10sp, contrast=WB_10sp.contrasts[,"Tpa_S_A"])
Tps_lrt_WB_10sp_S_A  <- glmLRT(fit_WB_10sp, contrast=WB_10sp.contrasts[,"Tps_S_A"])


RT_10sp.contrasts <- makeContrasts(
Tbi_S_A	=	Tbi.A - Tbi.S ,   ### +ve FC = higher in ASEX 
Tce_S_A	=	Tce.A - Tce.S ,
Tcm_S_A	=	Tcm.A - Tcm.S ,
Tpa_S_A	=	Tpa.A - Tpa.S ,
Tps_S_A	=	Tps.A - Tps.S ,
levels=mat_RT_10sp)

Tbi_lrt_RT_10sp_S_A  <- glmLRT(fit_RT_10sp, contrast=RT_10sp.contrasts[,"Tbi_S_A"])
Tce_lrt_RT_10sp_S_A  <- glmLRT(fit_RT_10sp, contrast=RT_10sp.contrasts[,"Tce_S_A"])
Tcm_lrt_RT_10sp_S_A  <- glmLRT(fit_RT_10sp, contrast=RT_10sp.contrasts[,"Tcm_S_A"])
Tpa_lrt_RT_10sp_S_A  <- glmLRT(fit_RT_10sp, contrast=RT_10sp.contrasts[,"Tpa_S_A"])
Tps_lrt_RT_10sp_S_A  <- glmLRT(fit_RT_10sp, contrast=RT_10sp.contrasts[,"Tps_S_A"])


LG_10sp.contrasts <- makeContrasts(
Tbi_S_A	=	Tbi.A - Tbi.S ,   ### +ve FC = higher in ASEX 
Tce_S_A	=	Tce.A - Tce.S ,
Tcm_S_A	=	Tcm.A - Tcm.S ,
Tpa_S_A	=	Tpa.A - Tpa.S ,
Tps_S_A	=	Tps.A - Tps.S ,
levels=mat_LG_10sp)

Tbi_lrt_LG_10sp_S_A  <- glmLRT(fit_LG_10sp, contrast=LG_10sp.contrasts[,"Tbi_S_A"])
Tce_lrt_LG_10sp_S_A  <- glmLRT(fit_LG_10sp, contrast=LG_10sp.contrasts[,"Tce_S_A"])
Tcm_lrt_LG_10sp_S_A  <- glmLRT(fit_LG_10sp, contrast=LG_10sp.contrasts[,"Tcm_S_A"])
Tpa_lrt_LG_10sp_S_A  <- glmLRT(fit_LG_10sp, contrast=LG_10sp.contrasts[,"Tpa_S_A"])
Tps_lrt_LG_10sp_S_A  <- glmLRT(fit_LG_10sp, contrast=LG_10sp.contrasts[,"Tps_S_A"])


#################################################################################################
### get DE_gene_table (samples)

TTT_Tbi_lrt_WB_10sp_S_A <- get_DE_genes(Tbi_lrt_WB_10sp_S_A, 0.05)
TTT_Tce_lrt_WB_10sp_S_A <- get_DE_genes(Tce_lrt_WB_10sp_S_A, 0.05)
TTT_Tcm_lrt_WB_10sp_S_A <- get_DE_genes(Tcm_lrt_WB_10sp_S_A, 0.05)
TTT_Tpa_lrt_WB_10sp_S_A <- get_DE_genes(Tpa_lrt_WB_10sp_S_A, 0.05)
TTT_Tps_lrt_WB_10sp_S_A <- get_DE_genes(Tps_lrt_WB_10sp_S_A, 0.05)

TTT_Tbi_lrt_RT_10sp_S_A <- get_DE_genes(Tbi_lrt_RT_10sp_S_A, 0.05)
TTT_Tce_lrt_RT_10sp_S_A <- get_DE_genes(Tce_lrt_RT_10sp_S_A, 0.05)
TTT_Tcm_lrt_RT_10sp_S_A <- get_DE_genes(Tcm_lrt_RT_10sp_S_A, 0.05)
TTT_Tpa_lrt_RT_10sp_S_A <- get_DE_genes(Tpa_lrt_RT_10sp_S_A, 0.05)
TTT_Tps_lrt_RT_10sp_S_A <- get_DE_genes(Tps_lrt_RT_10sp_S_A, 0.05)

TTT_Tbi_lrt_LG_10sp_S_A <- get_DE_genes(Tbi_lrt_LG_10sp_S_A, 0.05)
TTT_Tce_lrt_LG_10sp_S_A <- get_DE_genes(Tce_lrt_LG_10sp_S_A, 0.05)
TTT_Tcm_lrt_LG_10sp_S_A <- get_DE_genes(Tcm_lrt_LG_10sp_S_A, 0.05)
TTT_Tpa_lrt_LG_10sp_S_A <- get_DE_genes(Tpa_lrt_LG_10sp_S_A, 0.05)
TTT_Tps_lrt_LG_10sp_S_A <- get_DE_genes(Tps_lrt_LG_10sp_S_A, 0.05)

####### export full tables

write.table(TTT_Tbi_lrt_WB_10sp_S_A$table, "10sp_DE_output/TTT_Tbi_lrt_WB_10sp_S_A.csv", sep = ',', quote = FALSE, row.names = FALSE)
write.table(TTT_Tce_lrt_WB_10sp_S_A$table, "10sp_DE_output/TTT_Tce_lrt_WB_10sp_S_A.csv", sep = ',', quote = FALSE, row.names = FALSE)
write.table(TTT_Tcm_lrt_WB_10sp_S_A$table, "10sp_DE_output/TTT_Tcm_lrt_WB_10sp_S_A.csv", sep = ',', quote = FALSE, row.names = FALSE)
write.table(TTT_Tpa_lrt_WB_10sp_S_A$table, "10sp_DE_output/TTT_Tpa_lrt_WB_10sp_S_A.csv", sep = ',', quote = FALSE, row.names = FALSE)
write.table(TTT_Tps_lrt_WB_10sp_S_A$table, "10sp_DE_output/TTT_Tps_lrt_WB_10sp_S_A.csv", sep = ',', quote = FALSE, row.names = FALSE)
write.table(TTT_Tbi_lrt_RT_10sp_S_A$table, "10sp_DE_output/TTT_Tbi_lrt_RT_10sp_S_A.csv", sep = ',', quote = FALSE, row.names = FALSE)
write.table(TTT_Tce_lrt_RT_10sp_S_A$table, "10sp_DE_output/TTT_Tce_lrt_RT_10sp_S_A.csv", sep = ',', quote = FALSE, row.names = FALSE)
write.table(TTT_Tcm_lrt_RT_10sp_S_A$table, "10sp_DE_output/TTT_Tcm_lrt_RT_10sp_S_A.csv", sep = ',', quote = FALSE, row.names = FALSE)
write.table(TTT_Tpa_lrt_RT_10sp_S_A$table, "10sp_DE_output/TTT_Tpa_lrt_RT_10sp_S_A.csv", sep = ',', quote = FALSE, row.names = FALSE)
write.table(TTT_Tps_lrt_RT_10sp_S_A$table, "10sp_DE_output/TTT_Tps_lrt_RT_10sp_S_A.csv", sep = ',', quote = FALSE, row.names = FALSE)
write.table(TTT_Tbi_lrt_LG_10sp_S_A$table, "10sp_DE_output/TTT_Tbi_lrt_LG_10sp_S_A.csv", sep = ',', quote = FALSE, row.names = FALSE)
write.table(TTT_Tce_lrt_LG_10sp_S_A$table, "10sp_DE_output/TTT_Tce_lrt_LG_10sp_S_A.csv", sep = ',', quote = FALSE, row.names = FALSE)
write.table(TTT_Tcm_lrt_LG_10sp_S_A$table, "10sp_DE_output/TTT_Tcm_lrt_LG_10sp_S_A.csv", sep = ',', quote = FALSE, row.names = FALSE)
write.table(TTT_Tpa_lrt_LG_10sp_S_A$table, "10sp_DE_output/TTT_Tpa_lrt_LG_10sp_S_A.csv", sep = ',', quote = FALSE, row.names = FALSE)
write.table(TTT_Tps_lrt_LG_10sp_S_A$table, "10sp_DE_output/TTT_Tps_lrt_LG_10sp_S_A.csv", sep = ',', quote = FALSE, row.names = FALSE)

#######################################################################################################################

## N genes plot

N_sig_genes_10sp_S_A <- as.data.frame(c(
length(TTT_Tbi_lrt_WB_10sp_S_A$S_gene_list),
length(TTT_Tce_lrt_WB_10sp_S_A$S_gene_list),
length(TTT_Tcm_lrt_WB_10sp_S_A$S_gene_list),
length(TTT_Tpa_lrt_WB_10sp_S_A$S_gene_list),
length(TTT_Tps_lrt_WB_10sp_S_A$S_gene_list),

length(TTT_Tbi_lrt_RT_10sp_S_A$S_gene_list),
length(TTT_Tce_lrt_RT_10sp_S_A$S_gene_list),
length(TTT_Tcm_lrt_RT_10sp_S_A$S_gene_list),
length(TTT_Tpa_lrt_RT_10sp_S_A$S_gene_list),
length(TTT_Tps_lrt_RT_10sp_S_A$S_gene_list),

length(TTT_Tbi_lrt_LG_10sp_S_A$S_gene_list),
length(TTT_Tce_lrt_LG_10sp_S_A$S_gene_list),
length(TTT_Tcm_lrt_LG_10sp_S_A$S_gene_list),
length(TTT_Tpa_lrt_LG_10sp_S_A$S_gene_list),
length(TTT_Tps_lrt_LG_10sp_S_A$S_gene_list)
))

names(N_sig_genes_10sp_S_A) <- "Nsiggenes"
N_sig_genes_10sp_S_A$sp <- c("Tbi","Tce","Tcm","Tpa","Tps","Tbi","Tce","Tcm","Tpa","Tps","Tbi","Tce","Tcm","Tpa","Tps")
N_sig_genes_10sp_S_A$tiss <- c("WB","WB","WB","WB","WB","RT","RT","RT","RT","RT","LG","LG","LG","LG","LG")
N_sig_genes_10sp_S_A$group <- paste(N_sig_genes_10sp_S_A$sp, N_sig_genes_10sp_S_A$tiss, sep = "_")

N_sig_genes_10sp_S_A$group_ordered <-  ordered(N_sig_genes_10sp_S_A$group, levels=c(
"Tbi_WB","Tbi_RT","Tbi_LG",
"Tce_WB","Tce_RT","Tce_LG",
"Tps_WB","Tps_RT","Tps_LG",
"Tcm_WB","Tcm_RT","Tcm_LG",
"Tpa_WB","Tpa_RT","Tpa_LG"
))

p3 <- ggplot(N_sig_genes_10sp_S_A, aes(factor(group_ordered), Nsiggenes, fill = tiss )) + 
	geom_bar(stat="identity", position = "dodge") + 
	theme_bw() +
	xlab ("Species pair") + 
	ylab ("Number of DE genes, FDR < 0.05")  +
	scale_y_continuous(expand = c(0,0)) + 
	coord_cartesian(ylim=c(-0,650)) 
	
p3 + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
   scale_fill_manual(values=c("#56B4E9", "firebrick3", "black"))



###### filt by log fold change

FC_thesh = 1

N_sig_genes_10sp_S_A_logFC_1 <- as.data.frame(c(
sum(sqrt(TTT_Tbi_lrt_WB_10sp_S_A$S_logFC_list * TTT_Tbi_lrt_WB_10sp_S_A$S_logFC_list) > FC_thesh),
sum(sqrt(TTT_Tce_lrt_WB_10sp_S_A$S_logFC_list * TTT_Tce_lrt_WB_10sp_S_A$S_logFC_list) > FC_thesh),
sum(sqrt(TTT_Tcm_lrt_WB_10sp_S_A$S_logFC_list * TTT_Tcm_lrt_WB_10sp_S_A$S_logFC_list) > FC_thesh),
sum(sqrt(TTT_Tpa_lrt_WB_10sp_S_A$S_logFC_list * TTT_Tpa_lrt_WB_10sp_S_A$S_logFC_list) > FC_thesh),
sum(sqrt(TTT_Tps_lrt_WB_10sp_S_A$S_logFC_list * TTT_Tps_lrt_WB_10sp_S_A$S_logFC_list) > FC_thesh),

sum(sqrt(TTT_Tbi_lrt_RT_10sp_S_A$S_logFC_list * TTT_Tbi_lrt_RT_10sp_S_A$S_logFC_list) > FC_thesh),
sum(sqrt(TTT_Tce_lrt_RT_10sp_S_A$S_logFC_list * TTT_Tce_lrt_RT_10sp_S_A$S_logFC_list) > FC_thesh),
sum(sqrt(TTT_Tcm_lrt_RT_10sp_S_A$S_logFC_list * TTT_Tcm_lrt_RT_10sp_S_A$S_logFC_list) > FC_thesh),
sum(sqrt(TTT_Tpa_lrt_RT_10sp_S_A$S_logFC_list * TTT_Tpa_lrt_RT_10sp_S_A$S_logFC_list) > FC_thesh),
sum(sqrt(TTT_Tps_lrt_RT_10sp_S_A$S_logFC_list * TTT_Tps_lrt_RT_10sp_S_A$S_logFC_list) > FC_thesh),

sum(sqrt(TTT_Tbi_lrt_LG_10sp_S_A$S_logFC_list * TTT_Tbi_lrt_LG_10sp_S_A$S_logFC_list) > FC_thesh),
sum(sqrt(TTT_Tce_lrt_LG_10sp_S_A$S_logFC_list * TTT_Tce_lrt_LG_10sp_S_A$S_logFC_list) > FC_thesh),
sum(sqrt(TTT_Tcm_lrt_LG_10sp_S_A$S_logFC_list * TTT_Tcm_lrt_LG_10sp_S_A$S_logFC_list) > FC_thesh),
sum(sqrt(TTT_Tpa_lrt_LG_10sp_S_A$S_logFC_list * TTT_Tpa_lrt_LG_10sp_S_A$S_logFC_list) > FC_thesh),
sum(sqrt(TTT_Tps_lrt_LG_10sp_S_A$S_logFC_list * TTT_Tps_lrt_LG_10sp_S_A$S_logFC_list) > FC_thesh)
))


names(N_sig_genes_10sp_S_A_logFC_1) <- "Nsiggenes"
N_sig_genes_10sp_S_A_logFC_1$sp <- c("Tbi","Tce","Tcm","Tpa","Tps","Tbi","Tce","Tcm","Tpa","Tps","Tbi","Tce","Tcm","Tpa","Tps")
N_sig_genes_10sp_S_A_logFC_1$tiss <- c("WB","WB","WB","WB","WB","RT","RT","RT","RT","RT","LG","LG","LG","LG","LG")
N_sig_genes_10sp_S_A_logFC_1$group <- paste(N_sig_genes_10sp_S_A$sp, N_sig_genes_10sp_S_A$tiss, sep = "_")

N_sig_genes_10sp_S_A_logFC_1$group_ordered <-  ordered(N_sig_genes_10sp_S_A_logFC_1$group, levels=c(
"Tbi_WB","Tbi_RT","Tbi_LG",
"Tce_WB","Tce_RT","Tce_LG",
"Tps_WB","Tps_RT","Tps_LG",
"Tcm_WB","Tcm_RT","Tcm_LG",
"Tpa_WB","Tpa_RT","Tpa_LG"
))

p3a <- ggplot(N_sig_genes_10sp_S_A_logFC_1, aes(factor(group_ordered), Nsiggenes, fill = tiss )) + 
	geom_bar(stat="identity", position = "dodge") + 
	theme_bw() +
	xlab ("Species pair") + 
	ylab ("Number of DE genes, FDR < 0.05, FC > 2")  +
	scale_y_continuous(expand = c(0,0)) + 
	coord_cartesian(ylim=c(-0,300)) 
	
p3a + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
   scale_fill_manual(values=c("#56B4E9", "firebrick3", "black"))


FC_thesh = 2

N_sig_genes_10sp_S_A_logFC_2 <- as.data.frame(c(
sum(sqrt(TTT_Tbi_lrt_WB_10sp_S_A$S_logFC_list * TTT_Tbi_lrt_WB_10sp_S_A$S_logFC_list) > FC_thesh),
sum(sqrt(TTT_Tce_lrt_WB_10sp_S_A$S_logFC_list * TTT_Tce_lrt_WB_10sp_S_A$S_logFC_list) > FC_thesh),
sum(sqrt(TTT_Tcm_lrt_WB_10sp_S_A$S_logFC_list * TTT_Tcm_lrt_WB_10sp_S_A$S_logFC_list) > FC_thesh),
sum(sqrt(TTT_Tpa_lrt_WB_10sp_S_A$S_logFC_list * TTT_Tpa_lrt_WB_10sp_S_A$S_logFC_list) > FC_thesh),
sum(sqrt(TTT_Tps_lrt_WB_10sp_S_A$S_logFC_list * TTT_Tps_lrt_WB_10sp_S_A$S_logFC_list) > FC_thesh),

sum(sqrt(TTT_Tbi_lrt_RT_10sp_S_A$S_logFC_list * TTT_Tbi_lrt_RT_10sp_S_A$S_logFC_list) > FC_thesh),
sum(sqrt(TTT_Tce_lrt_RT_10sp_S_A$S_logFC_list * TTT_Tce_lrt_RT_10sp_S_A$S_logFC_list) > FC_thesh),
sum(sqrt(TTT_Tcm_lrt_RT_10sp_S_A$S_logFC_list * TTT_Tcm_lrt_RT_10sp_S_A$S_logFC_list) > FC_thesh),
sum(sqrt(TTT_Tpa_lrt_RT_10sp_S_A$S_logFC_list * TTT_Tpa_lrt_RT_10sp_S_A$S_logFC_list) > FC_thesh),
sum(sqrt(TTT_Tps_lrt_RT_10sp_S_A$S_logFC_list * TTT_Tps_lrt_RT_10sp_S_A$S_logFC_list) > FC_thesh),

sum(sqrt(TTT_Tbi_lrt_LG_10sp_S_A$S_logFC_list * TTT_Tbi_lrt_LG_10sp_S_A$S_logFC_list) > FC_thesh),
sum(sqrt(TTT_Tce_lrt_LG_10sp_S_A$S_logFC_list * TTT_Tce_lrt_LG_10sp_S_A$S_logFC_list) > FC_thesh),
sum(sqrt(TTT_Tcm_lrt_LG_10sp_S_A$S_logFC_list * TTT_Tcm_lrt_LG_10sp_S_A$S_logFC_list) > FC_thesh),
sum(sqrt(TTT_Tpa_lrt_LG_10sp_S_A$S_logFC_list * TTT_Tpa_lrt_LG_10sp_S_A$S_logFC_list) > FC_thesh),
sum(sqrt(TTT_Tps_lrt_LG_10sp_S_A$S_logFC_list * TTT_Tps_lrt_LG_10sp_S_A$S_logFC_list) > FC_thesh)
))


names(N_sig_genes_10sp_S_A_logFC_2) <- "Nsiggenes"
N_sig_genes_10sp_S_A_logFC_2$sp <- c("Tbi","Tce","Tcm","Tpa","Tps","Tbi","Tce","Tcm","Tpa","Tps","Tbi","Tce","Tcm","Tpa","Tps")
N_sig_genes_10sp_S_A_logFC_2$tiss <- c("WB","WB","WB","WB","WB","RT","RT","RT","RT","RT","LG","LG","LG","LG","LG")
N_sig_genes_10sp_S_A_logFC_2$group <- paste(N_sig_genes_10sp_S_A$sp, N_sig_genes_10sp_S_A$tiss, sep = "_")

N_sig_genes_10sp_S_A_logFC_2$group_ordered <-  ordered(N_sig_genes_10sp_S_A_logFC_1$group, levels=c(
"Tbi_WB","Tbi_RT","Tbi_LG",
"Tce_WB","Tce_RT","Tce_LG",
"Tps_WB","Tps_RT","Tps_LG",
"Tcm_WB","Tcm_RT","Tcm_LG",
"Tpa_WB","Tpa_RT","Tpa_LG"
))

p3b <- ggplot(N_sig_genes_10sp_S_A_logFC_2, aes(factor(group_ordered), Nsiggenes, fill = tiss )) + 
	geom_bar(stat="identity", position = "dodge") + 
	theme_bw() +
	xlab ("Species pair") + 
	ylab ("Number of DE genes, FDR < 0.05, FC > 4")  +
	scale_y_continuous(expand = c(0,0)) + 
	coord_cartesian(ylim=c(-0,80)) 
	
p3b + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
   scale_fill_manual(values=c("#56B4E9", "firebrick3", "black"))


##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### lets venn

##### diff sp code


five_sp_DE_venn <- function(Tbi,Tce,Tcm,Tpa,Tps,title){
	venny.plot <- venn.diagram(
	list("Tbi" = Tbi, "Tce" = Tce, "Tcm" = Tcm, "Tpa" = Tpa, "Tps" = Tps ), filename = NULL,
                            fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
                            cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
                            margin = 0.6, cat.dist = 0.23, main = title, main.pos = c(0.5,0.8), main.cex = 2, main.fontface = "bold", cat.cex = 2)
	return(venny.plot)
}




#### diff sp samples

WB_10sp_SvsA_N_sig_genes <- five_sp_DE_venn(
TTT_Tbi_lrt_WB_10sp_S_A$S_gene_list,
TTT_Tce_lrt_WB_10sp_S_A$S_gene_list,
TTT_Tcm_lrt_WB_10sp_S_A$S_gene_list,
TTT_Tpa_lrt_WB_10sp_S_A$S_gene_list,
TTT_Tps_lrt_WB_10sp_S_A$S_gene_list, 
"10 sp orths, whole-body, FDR = 0.05"
)

RT_10sp_SvsA_N_sig_genes <- five_sp_DE_venn(
TTT_Tbi_lrt_RT_10sp_S_A$S_gene_list,
TTT_Tce_lrt_RT_10sp_S_A$S_gene_list,
TTT_Tcm_lrt_RT_10sp_S_A$S_gene_list,
TTT_Tpa_lrt_RT_10sp_S_A$S_gene_list,
TTT_Tps_lrt_RT_10sp_S_A$S_gene_list,
"10 sp orths, rep tract, FDR = 0.05"
)


LG_10sp_SvsA_N_sig_genes <- five_sp_DE_venn(
TTT_Tbi_lrt_LG_10sp_S_A$S_gene_list,
TTT_Tce_lrt_LG_10sp_S_A$S_gene_list,
TTT_Tcm_lrt_LG_10sp_S_A$S_gene_list,
TTT_Tpa_lrt_LG_10sp_S_A$S_gene_list,
TTT_Tps_lrt_LG_10sp_S_A$S_gene_list, 
"10 sp orths, Legs, FDR = 0.05"
)

grid.arrange(gTree(children=WB_10sp_SvsA_N_sig_genes), gTree(children=RT_10sp_SvsA_N_sig_genes), gTree(children=LG_10sp_SvsA_N_sig_genes))

################################################################################################################################
#### calc prob of overlaps ### note these produce uncorrected pvals

WB_10sp_SvsA_N_sig_genes_list <- list(
TTT_Tps_lrt_WB_10sp_S_A$S_gene_list,
TTT_Tcm_lrt_WB_10sp_S_A$S_gene_list,
TTT_Tce_lrt_WB_10sp_S_A$S_gene_list,
TTT_Tpa_lrt_WB_10sp_S_A$S_gene_list,
TTT_Tbi_lrt_WB_10sp_S_A$S_gene_list
) 

str(WB_10sp_SvsA_N_sig_genes_list)
WB_10sp_SvsA_N_total_genes = length(t(y_WB_10sp$genes))
WB_10sp_SvsA_N_total_genes


sup_t_WB_10sp_SvsA = supertest(WB_10sp_SvsA_N_sig_genes_list, n=WB_10sp_SvsA_N_total_genes)

#######
RT_10sp_SvsA_N_sig_genes_list <- list(
TTT_Tps_lrt_RT_10sp_S_A$S_gene_list,
TTT_Tcm_lrt_RT_10sp_S_A$S_gene_list,
TTT_Tce_lrt_RT_10sp_S_A$S_gene_list,
TTT_Tpa_lrt_RT_10sp_S_A$S_gene_list,
TTT_Tbi_lrt_RT_10sp_S_A$S_gene_list) 

str(RT_10sp_SvsA_N_sig_genes_list)
RT_10sp_SvsA_N_total_genes = length(t(y_RT_10sp$genes))
RT_10sp_SvsA_N_total_genes

sup_t_RT_10sp_SvsA = supertest(RT_10sp_SvsA_N_sig_genes_list, n=RT_10sp_SvsA_N_total_genes)


#######
LG_10sp_SvsA_N_sig_genes_list <- list(
TTT_Tps_lrt_LG_10sp_S_A$S_gene_list,
TTT_Tcm_lrt_LG_10sp_S_A$S_gene_list,
TTT_Tce_lrt_LG_10sp_S_A$S_gene_list,
TTT_Tpa_lrt_LG_10sp_S_A$S_gene_list,
TTT_Tbi_lrt_LG_10sp_S_A$S_gene_list) 

str(LG_10sp_SvsA_N_sig_genes_list)
LG_10sp_SvsA_N_total_genes = length(t(y_LG_10sp$genes))
LG_10sp_SvsA_N_total_genes

sup_t_LG_10sp_SvsA = supertest(LG_10sp_SvsA_N_sig_genes_list, n=LG_10sp_SvsA_N_total_genes)

### export tables (drawing here is a pain)

write.csv(summary(sup_t_WB_10sp_SvsA)$Table, file="10sp_DE_output/sup_t_WB_10sp_SvsA_summary.table.csv", row.names=FALSE)
write.csv(summary(sup_t_RT_10sp_SvsA)$Table, file="10sp_DE_output/sup_t_RT_10sp_SvsA_summary.table.csv", row.names=FALSE)
write.csv(summary(sup_t_LG_10sp_SvsA)$Table, file="10sp_DE_output/sup_t_LG_10sp_SvsA_summary.table.csv", row.names=FALSE)


#############################################################################################################
## diff tiss overlap


tiss_WBRTLG_venn <- function(WB,RT,LG,title){
	venny.plot <- venn.diagram(
	list("WB" = WB, "RT" = RT, "LG" = LG), filename = NULL,
                            fill    = c("dodgerblue", "goldenrod1", "seagreen3"),
                            cat.col = c("dodgerblue", "goldenrod1", "seagreen3"),
                            margin = 0.6, cat.dist = 0.05, main = title, main.pos = c(0.5,0.75), main.cex = 2, main.fontface = "bold", cat.cex = 2, cex = 2)
	return(venny.plot)
}

Tbi_WBRTLG_10sp_SvsA_N_sig_genes <- tiss_WBRTLG_venn(
TTT_Tbi_lrt_WB_10sp_S_A$S_gene_list,
TTT_Tbi_lrt_RT_10sp_S_A$S_gene_list,
TTT_Tbi_lrt_LG_10sp_S_A$S_gene_list,
"10 sp orths, Tbi, FDR = 0.05"
)

Tce_WBRTLG_10sp_SvsA_N_sig_genes <- tiss_WBRTLG_venn(
TTT_Tce_lrt_WB_10sp_S_A$S_gene_list,
TTT_Tce_lrt_RT_10sp_S_A$S_gene_list,
TTT_Tce_lrt_LG_10sp_S_A$S_gene_list,
"10 sp orths, Tce, FDR = 0.05"
)

Tcm_WBRTLG_10sp_SvsA_N_sig_genes <- tiss_WBRTLG_venn(
TTT_Tcm_lrt_WB_10sp_S_A$S_gene_list,
TTT_Tcm_lrt_RT_10sp_S_A$S_gene_list,
TTT_Tcm_lrt_LG_10sp_S_A$S_gene_list,
"10 sp orths, Tcm, FDR = 0.05"
)

Tpa_WBRTLG_10sp_SvsA_N_sig_genes <- tiss_WBRTLG_venn(
TTT_Tpa_lrt_WB_10sp_S_A$S_gene_list,
TTT_Tpa_lrt_RT_10sp_S_A$S_gene_list,
TTT_Tpa_lrt_LG_10sp_S_A$S_gene_list,
"10 sp orths, Tpa, FDR = 0.05"
)

Tps_WBRTLG_10sp_SvsA_N_sig_genes <- tiss_WBRTLG_venn(
TTT_Tps_lrt_WB_10sp_S_A$S_gene_list,
TTT_Tps_lrt_RT_10sp_S_A$S_gene_list,
TTT_Tps_lrt_LG_10sp_S_A$S_gene_list,
"10 sp orths, Tps, FDR = 0.05"
)

grid.arrange(gTree(children=Tbi_WBRTLG_10sp_SvsA_N_sig_genes),gTree(children=Tce_WBRTLG_10sp_SvsA_N_sig_genes),
                gTree(children=Tcm_WBRTLG_10sp_SvsA_N_sig_genes),gTree(children=Tpa_WBRTLG_10sp_SvsA_N_sig_genes),
                gTree(children=Tps_WBRTLG_10sp_SvsA_N_sig_genes))


### superexact 

## Tbi
Tbi_WBRTLG_tiss_sig_st_list <- list(
TTT_Tbi_lrt_WB_10sp_S_A$S_gene_list,
TTT_Tbi_lrt_RT_10sp_S_A$S_gene_list,
TTT_Tbi_lrt_LG_10sp_S_A$S_gene_list)

## get lens
### lens differ (as differnt genes exculded in diff tissues) 
### to make more conservative I will set the background length to the inntersect size of all tissues

Tbi_WBRTLG_tiss_lens_all <- list(
WB = TTT_Tbi_lrt_WB_10sp_S_A$table$genes,
RT = TTT_Tbi_lrt_RT_10sp_S_A$table$genes,
LG = TTT_Tbi_lrt_LG_10sp_S_A$table$genes)

Tbi_WBRTLG_tiss_lens_intersect = length(Intersect(Tbi_WBRTLG_tiss_lens_all ))
Tbi_WBRTLG_tiss_lens_min_st = supertest(Tbi_WBRTLG_tiss_sig_st_list , n=Tbi_WBRTLG_tiss_lens_intersect)
summary(Tbi_WBRTLG_tiss_lens_min_st)

write.csv(summary(Tbi_WBRTLG_tiss_lens_min_st)$Table, file="10sp_DE_output/sup_Tbi_WBRTLG_tiss_lens_min_st_summary.table.csv", row.names=FALSE)

## Tce

Tce_WBRTLG_tiss_sig_st_list <- list(
TTT_Tce_lrt_WB_10sp_S_A$S_gene_list,
TTT_Tce_lrt_RT_10sp_S_A$S_gene_list,
TTT_Tce_lrt_LG_10sp_S_A$S_gene_list)

## get lens
### lens differ (as differnt genes exculded in diff tissues) 
### to make more conservative I will set the background length to the inntersect size of all tissues

Tce_WBRTLG_tiss_lens_all <- list(
WB = TTT_Tce_lrt_WB_10sp_S_A$table$genes,
RT = TTT_Tce_lrt_RT_10sp_S_A$table$genes,
LG = TTT_Tce_lrt_LG_10sp_S_A$table$genes)

Tce_WBRTLG_tiss_lens_intersect = length(Intersect(Tce_WBRTLG_tiss_lens_all ))

### superexact 
Tce_WBRTLG_tiss_lens_min_st = supertest(Tce_WBRTLG_tiss_sig_st_list , n=Tce_WBRTLG_tiss_lens_intersect)
summary(Tce_WBRTLG_tiss_lens_min_st)

write.csv(summary(Tce_WBRTLG_tiss_lens_min_st)$Table, file="10sp_DE_output/sup_Tce_WBRTLG_tiss_lens_min_st_summary.table.csv", row.names=FALSE)

## Tcm

Tcm_WBRTLG_tiss_sig_st_list <- list(
TTT_Tcm_lrt_WB_10sp_S_A$S_gene_list,
TTT_Tcm_lrt_RT_10sp_S_A$S_gene_list,
TTT_Tcm_lrt_LG_10sp_S_A$S_gene_list)

## get lens
### lens differ (as differnt genes exculded in diff tissues) 
### to make more conservative I will set the background length to the inntersect size of all tissues

Tcm_WBRTLG_tiss_lens_all <- list(
WB = TTT_Tcm_lrt_WB_10sp_S_A$table$genes,
RT = TTT_Tcm_lrt_RT_10sp_S_A$table$genes,
LG = TTT_Tcm_lrt_LG_10sp_S_A$table$genes)

Tcm_WBRTLG_tiss_lens_intersect = length(Intersect(Tcm_WBRTLG_tiss_lens_all ))

### superexact 
Tcm_WBRTLG_tiss_lens_min_st = supertest(Tcm_WBRTLG_tiss_sig_st_list , n=Tcm_WBRTLG_tiss_lens_intersect)
summary(Tcm_WBRTLG_tiss_lens_min_st)

write.csv(summary(Tcm_WBRTLG_tiss_lens_min_st)$Table, file="10sp_DE_output/sup_Tcm_WBRTLG_tiss_lens_min_st_summary.table.csv", row.names=FALSE)


## Tpa
Tpa_WBRTLG_tiss_sig_st_list <- list(
TTT_Tpa_lrt_WB_10sp_S_A$S_gene_list,
TTT_Tpa_lrt_RT_10sp_S_A$S_gene_list,
TTT_Tpa_lrt_LG_10sp_S_A$S_gene_list)

## get lens
### lens differ (as differnt genes exculded in diff tissues) 
### to make more conservative I will set the background length to the inntersect size of all tissues

Tpa_WBRTLG_tiss_lens_all <- list(
WB = TTT_Tpa_lrt_WB_10sp_S_A$table$genes,
RT = TTT_Tpa_lrt_RT_10sp_S_A$table$genes,
LG = TTT_Tpa_lrt_LG_10sp_S_A$table$genes)

Tpa_WBRTLG_tiss_lens_intersect = length(Intersect(Tpa_WBRTLG_tiss_lens_all ))

### superexact 
Tpa_WBRTLG_tiss_lens_min_st = supertest(Tpa_WBRTLG_tiss_sig_st_list , n=Tpa_WBRTLG_tiss_lens_intersect)
summary(Tpa_WBRTLG_tiss_lens_min_st)

write.csv(summary(Tpa_WBRTLG_tiss_lens_min_st)$Table, file="10sp_DE_output/sup_Tpa_WBRTLG_tiss_lens_min_st_summary.table.csv", row.names=FALSE)

## Tps

Tps_WBRTLG_tiss_sig_st_list <- list(
TTT_Tps_lrt_WB_10sp_S_A$S_gene_list,
TTT_Tps_lrt_RT_10sp_S_A$S_gene_list,
TTT_Tps_lrt_LG_10sp_S_A$S_gene_list)

## get lens
### lens differ (as differnt genes exculded in diff tissues) 
### to make more conservative I will set the background length to the inntersect size of all tissues

Tps_WBRTLG_tiss_lens_all <- list(
WB = TTT_Tps_lrt_WB_10sp_S_A$table$genes,
RT = TTT_Tps_lrt_RT_10sp_S_A$table$genes,
LG = TTT_Tps_lrt_LG_10sp_S_A$table$genes)

Tps_WBRTLG_tiss_lens_intersect = length(Intersect(Tps_WBRTLG_tiss_lens_all ))

### superexact 
Tps_WBRTLG_tiss_lens_min_st = supertest(Tps_WBRTLG_tiss_sig_st_list , n=Tps_WBRTLG_tiss_lens_intersect)
summary(Tps_WBRTLG_tiss_lens_min_st)

write.csv(summary(Tps_WBRTLG_tiss_lens_min_st)$Table, file="10sp_DE_output/sup_Tps_WBRTLG_tiss_lens_min_st_summary.table.csv", row.names=FALSE)



###################################################################################################################################################################################
###### are DE genes also interaction genes?

### unique se DE genes (from contasts)
# wb
all_DE_genes_WB_10sp <- unique(c(
TTT_Tbi_lrt_WB_10sp_S_A$S_gene_list,
TTT_Tce_lrt_WB_10sp_S_A$S_gene_list,
TTT_Tcm_lrt_WB_10sp_S_A$S_gene_list,
TTT_Tpa_lrt_WB_10sp_S_A$S_gene_list,
TTT_Tps_lrt_WB_10sp_S_A$S_gene_list
))

# rt
all_DE_genes_RT_10sp <- unique(c(
TTT_Tbi_lrt_RT_10sp_S_A$S_gene_list,
TTT_Tce_lrt_RT_10sp_S_A$S_gene_list,
TTT_Tcm_lrt_RT_10sp_S_A$S_gene_list,
TTT_Tpa_lrt_RT_10sp_S_A$S_gene_list,
TTT_Tps_lrt_RT_10sp_S_A$S_gene_list
))


#lg
all_DE_genes_LG_10sp <- unique(c(
TTT_Tbi_lrt_LG_10sp_S_A$S_gene_list,
TTT_Tce_lrt_LG_10sp_S_A$S_gene_list,
TTT_Tcm_lrt_LG_10sp_S_A$S_gene_list,
TTT_Tpa_lrt_LG_10sp_S_A$S_gene_list,
TTT_Tps_lrt_LG_10sp_S_A$S_gene_list
))



### genes showing sig OVERALL S vs A effect

WB_10sp_S_A_overall_eff_genes <- TTT_WB_10sp_S_A_propeff$S_gene_list
RT_10sp_S_A_overall_eff_genes <- TTT_RT_10sp_S_A_propeff$S_gene_list
LG_10sp_S_A_overall_eff_genes <- TTT_LG_10sp_S_A_propeff$S_gene_list

### genes with sig interaction

WB_10sp_S_A_INTER_genes <- TTT_WB_10sp_S_A_INTER$S_gene_list
RT_10sp_S_A_INTER_genes <- TTT_RT_10sp_S_A_INTER$S_gene_list
LG_10sp_S_A_INTER_genes <- TTT_LG_10sp_S_A_INTER$S_gene_list


#### venns (code)
### Interaction in red, DE in grey

inter_DE_venn <- function(DE_vect,inter_vect,title){
	if(length(DE_vect) < length(inter_vect)){
		invert_deg = 180
		print("inverting")	
	} else {
		invert_deg = 0
		print("not inverting")	
	}
	
	venny.plot <- venn.diagram(
	list("DE" = DE_vect, "Interaction" = inter_vect), filename = NULL,
                            cat.col = c( "black",   "red"),
                            fill=c("black",   "red"), margin = 0.2, main = title, main.pos = c(0.5,0.8), cat.dist = 0.03, main.cex = 2, main.fontface = "bold", rotation.degree = invert_deg, cat.cex = 0, cex = 1.5)
	return(venny.plot)
}

#### to contrast S vs A
WB_10sp_inter_DE_venn <- inter_DE_venn(all_DE_genes_WB_10sp, WB_10sp_S_A_INTER_genes,"Whole-Body")
RT_10sp_inter_DE_venn <- inter_DE_venn(all_DE_genes_RT_10sp, RT_10sp_S_A_INTER_genes,"rep tract")
LG_10sp_inter_DE_venn <- inter_DE_venn(all_DE_genes_LG_10sp, LG_10sp_S_A_INTER_genes,"Legs")


grid.arrange(gTree(children=WB_10sp_inter_DE_venn), gTree(children=RT_10sp_inter_DE_venn) ,gTree(children=LG_10sp_inter_DE_venn  ),
top = "Black = sig [sex vs asex] effect in any contrast. Red = sig [lineage] * [sex vs asex] ")


#### to overall eff S vs A

WB_10sp_inter_DE_venn_overall_effSA <- inter_DE_venn(WB_10sp_S_A_overall_eff_genes, WB_10sp_S_A_INTER_genes,"Whole-Body")
RT_10sp_inter_DE_venn_overall_effSA <- inter_DE_venn(RT_10sp_S_A_overall_eff_genes , RT_10sp_S_A_INTER_genes,"Gonads")
LG_10sp_inter_DE_venn_overall_effSA <- inter_DE_venn(LG_10sp_S_A_overall_eff_genes , LG_10sp_S_A_INTER_genes,"Legs")

grid.arrange(gTree(children=WB_10sp_inter_DE_venn_overall_effSA), gTree(children=RT_10sp_inter_DE_venn_overall_effSA) ,gTree(children=LG_10sp_inter_DE_venn_overall_effSA ), 
top = "Black = Sig overall [sex vs asex] effect. Red = sig [lineage] * [sex vs asex] ")



#########################################################################################################
## extract overall effect genes vs no effect


Intersect <- function (x) {  
  # Multiple set version of intersect
  # x is a list
  if (length(x) == 1) {
    unlist(x)
  } else if (length(x) == 2) {
    intersect(x[[1]], x[[2]])
  } else if (length(x) > 2){
    intersect(x[[1]], Intersect(x[-1]))
  }
}

Union <- function (x) {  
  # Multiple set version of union
  # x is a list
  if (length(x) == 1) {
    unlist(x)
  } else if (length(x) == 2) {
    union(x[[1]], x[[2]])
  } else if (length(x) > 2) {
    union(x[[1]], Union(x[-1]))
  }
}

Setdiff <- function (x, y) {
  # Remove the union of the y's from the common x's. 
  # x and y are lists of characters.
  xx <- Intersect(x)
  yy <- Union(y)
  setdiff(xx, yy)
}



WB_10sp_inter_DE_venn_overall_effSA_list <- list(WB_ove_eff = WB_10sp_S_A_overall_eff_genes, WB_inter = WB_10sp_S_A_INTER_genes)
RT_10sp_inter_DE_venn_overall_effSA_list <- list(RT_ove_eff = RT_10sp_S_A_overall_eff_genes, RT_inter = RT_10sp_S_A_INTER_genes)
LG_10sp_inter_DE_venn_overall_effSA_list <- list(LG_ove_eff = LG_10sp_S_A_overall_eff_genes, LG_inter = LG_10sp_S_A_INTER_genes)

WB_ove_eff_genes <- Setdiff(WB_10sp_inter_DE_venn_overall_effSA_list["WB_ove_eff"], WB_10sp_inter_DE_venn_overall_effSA_list["WB_inter"])
RT_ove_eff_genes <- Setdiff(RT_10sp_inter_DE_venn_overall_effSA_list["RT_ove_eff"], RT_10sp_inter_DE_venn_overall_effSA_list["RT_inter"])
LG_ove_eff_genes <- Setdiff(LG_10sp_inter_DE_venn_overall_effSA_list["LG_ove_eff"], LG_10sp_inter_DE_venn_overall_effSA_list["LG_inter"])
length(WB_ove_eff_genes) 
length(RT_ove_eff_genes) 
length(LG_ove_eff_genes)


write.table(WB_ove_eff_genes, "10sp_DE_output/WB_convergent_genes.csv", sep = ',', quote = FALSE, row.names = FALSE)
write.table(RT_ove_eff_genes, "10sp_DE_output/RT_convergent_genes.csv", sep = ',', quote = FALSE, row.names = FALSE)
write.table(LG_ove_eff_genes, "10sp_DE_output/LG_convergent_genes.csv", sep = ',', quote = FALSE, row.names = FALSE)


######### how much overlap between are the convergent

conv_overlap <- tiss_WBRTLG_venn(
WB_ove_eff_genes,
RT_ove_eff_genes,
LG_ove_eff_genes,
"convergent overlap, FDR = 0.05"
)

grid.arrange(gTree(children=conv_overlap))


## super exact test

overall_eff_WBRTLG_tiss_sig_st_list <- list(
WB = WB_ove_eff_genes,
RT = RT_ove_eff_genes,
LG = LG_ove_eff_genes)

## get lens
### lens differ (as differnt genes exculded in diff tissues) 
### to make more conservative I will set the background length to the inntersect size of all tissues

overall_eff_WBRTLG_tiss_lens_all <- list(
WB = TTT_WB_10sp_S_A_INTER$table$genes,
RT = TTT_RT_10sp_S_A_INTER$table$genes,
LG = TTT_LG_10sp_S_A_INTER$table$genes)

overall_eff_WBRTLG_tiss_lens_intersect = length(Intersect(overall_eff_WBRTLG_tiss_lens_all))

### superexact 
overall_eff_WBRTLG_tiss_st = supertest(overall_eff_WBRTLG_tiss_sig_st_list, n=overall_eff_WBRTLG_tiss_lens_intersect)


write.csv(summary(overall_eff_WBRTLG_tiss_st)$Table, file="10sp_DE_output/sup_conv_WBRTLG_tiss_st_summary.table.csv", row.names=FALSE)


########################################################################################################################################################################
####### output session info
print (sessionInfo())
writeLines(capture.output(sessionInfo()), "10sp_DE_output/10sp_edgeR.R_sessionInfo.txt")











