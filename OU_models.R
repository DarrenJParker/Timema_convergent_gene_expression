## OU_models.R 
## Note most of this is based on the tutorial from Sam Price and Roi Holzman in UC Davis (http://www.eve.ucdavis.edu/~wainwrightlab/Roi/Site/Teaching.html)

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#### aborts if no command line args provided
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}


########################################################################################################################
# libs

library(ape)
library(ouch)

# print (sessionInfo())
# R version 3.4.1 (2017-06-30)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: macOS High Sierra 10.13.6

# Matrix products: default
# BLAS: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRblas.0.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRlapack.dylib

# locale:
# [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     

# other attached packages:
# [1] ouch_2.11-1   subplex_1.5-4 ape_5.2      

# loaded via a namespace (and not attached):
# [1] compiler_3.4.1  parallel_3.4.1  Rcpp_0.12.13    nlme_3.1-131    grid_3.4.1      lattice_0.20-35


########################################################################################################################
# data

gene_file_name <- args[1]
output_dir_name <- args[2]

dir.create(output_dir_name)

mytimtree <- read.nexus("Data/for_OU_models/Tim_tree.nex")
mytimdata <- read.table(gene_file_name)
## mytimdata <- read.table("./WB_GE_for_OU/OG-118.txt")

curr_gene_name <- colnames(mytimdata)[1]

setwd(output_dir_name)


########################################################################################################################
# Get tree into OUCH format using the ape2ouch function in ouch

myoutimtree<-ape2ouch(mytimtree)
#plot(myoutimtree)

## Now we want to get our data into a format that ouch can use.  

myoutimdataframe <- as(myoutimtree,"data.frame")
mytimdata$labels <- row.names(mytimdata)
myoutimdata <- merge(myoutimdataframe, mytimdata,by="labels",all=T)
row.names(myoutimdata) <- myoutimdata$nodes
newoutimtree <- ouchtree(nodes= myoutimdata$nodes, ancestors=myoutimdata$ancestors, times=myoutimdata$times, labels=myoutimdata$labels)

colnames(myoutimdata) <- c("labels", "nodes", "ancestors","times", "focal_gene")
myoutimdata$focal_gene <- as.numeric(myoutimdata$focal_gene)


########################################################################################################################
# First run a Brownian model with the brown function 
# This is model for what we would expect under drift

brown_focal_gene <- brown(data=myoutimdata["focal_gene"], tree=newoutimtree)
#summary(brown_focal_gene)


#######################################################################################################################
# Next fit the OU model 
# OU models are the brownian model with a parameter to specify optima to test if the observed values are more likely due to drift or selection to these optima

## selction model in this case is that Asex species have a different optima than sexuals.
myoutimdata$sel <-as.factor(
ifelse(myoutimdata$labels == "Tte", "asex", 
ifelse(myoutimdata$labels == "Tms", "asex", 
ifelse(myoutimdata$labels == "Tsi", "asex", 
ifelse(myoutimdata$labels == "Tge", "asex", 
ifelse(myoutimdata$labels == "Tdi", "asex", "sex"))))))

## fit
sel_focal_gene <- hansen(data=myoutimdata["focal_gene"],tree=newoutimtree, regimes=myoutimdata["sel"], sqrt.alpha=1, sigma=1, fit=TRUE)

#plot(sel_focal_gene) ## to check the optima are correct

#### output

oumodelcompare_tim_r <- as.data.frame(cbind(
summary(brown_focal_gene)$loglik, summary(brown_focal_gene)$dof, summary(sel_focal_gene)$loglik, summary(sel_focal_gene)$dof))

colnames(oumodelcompare_tim_r)<-c("Brown_loglik","Brown_df","OU_loglik","OU_df") 

### compare the models with an LRT

oumodelcompare_tim_r$LRT_OU_brown_p <- pchisq(2 * (oumodelcompare_tim_r$OU_loglik - oumodelcompare_tim_r$Brown_loglik), df=2, lower.tail=F)

oumodelcompare_tim_r_out <- cbind(curr_gene_name, oumodelcompare_tim_r)

#write.csv(oumodelcompare_tim_r_out, file=paste(curr_gene_name, "OUmodelcompare.csv", sep = "_"), quote = FALSE, row.names = FALSE)


#######################################################################################################################
# Next fit the OU model where the asex optima can be different to each other
# OU models are the brownian model with a parameter to specify optima to test if the observed values are more likely due to drift or selection to these optima

## selection model in this case is that Asex species have a different optima than sexuals and to each other.
myoutimdata$sel2 <-as.factor(
ifelse(myoutimdata$labels == "Tte", "asex1", 
ifelse(myoutimdata$labels == "Tms", "asex2", 
ifelse(myoutimdata$labels == "Tsi", "asex3", 
ifelse(myoutimdata$labels == "Tge", "asex4", 
ifelse(myoutimdata$labels == "Tdi", "asex5", "sex"))))))

## fit
sel_focal_gene2 <- hansen(data=myoutimdata["focal_gene"],tree=newoutimtree, regimes=myoutimdata["sel2"], sqrt.alpha=1, sigma=1, fit=TRUE)

#plot(sel_focal_gene2) ## to check the optima are correct

#### output

oumodelcompare_tim_r2 <- as.data.frame(cbind(oumodelcompare_tim_r,
summary(sel_focal_gene2)$loglik, summary(sel_focal_gene2)$dof))

colnames(oumodelcompare_tim_r2) <-c("Brown_loglik","Brown_df","OU_2opt_loglik","OU_2opt_df","LRT_OU_2opt_brown_p", "OU_multiopt_loglik","OU_multiopt_sel2_df") 

### compare the sel1 model with sel2 model with an LRT


oumodelcompare_tim_r2$LRT_OU_2opt_OU_multiopt_p <- pchisq(2 * (oumodelcompare_tim_r2$OU_multiopt_loglik - oumodelcompare_tim_r2$OU_2opt_loglik), df=4, lower.tail=F)


oumodelcompare_tim_r2_out <- cbind(curr_gene_name, oumodelcompare_tim_r2)

write.csv(oumodelcompare_tim_r2_out, file=paste(curr_gene_name, "OUmodelcompare.csv", sep = "_"), quote = FALSE, row.names = FALSE)


















