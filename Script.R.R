
# Importing packages required for all analysis

library(NOISeq)
library(DESeq2)
library(PoissonSeq)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(ggrepel)
library(clusterProfiler)
library(DEGreport)
library(org.Hs.eg.db)
library(DOSE)
library(pathview)
library(tidyverse)
library(EnsDb.Hsapiens.v86)
library(AnnotationHub)
library(ensembldb)
library(dplyr)
library(PROPER)

# Preparing the GEO data for analysis - Real Data
GSE40562 = getGEO("GSE40562")[[1]]
GSE98582 = getGEO("GSE98582")[[1]]
GSE40562_eset = exprs(GSE40562)
GSE98582_eset = exprs(GSE98582)

# Making the normalization methods from packages and coding by self
normalization_methods = list('DESeq2' = NA,
                             'PoissonSeq' = NA,
                             'TMM' = tmm, #edgeR and QuasiSeq
                             'RPKM' = rpkm,
                             'UQ' = uqua)

normalization_methods$DESeq2 = function(expression_set){
    temp = expression_set
    psuedo_ref_sample = col_multiply(expression_set)
    normalization_factor_matrix = sweep(expression_set,FUN="/",MARGIN=1,STATS=psuedo_ref_sample)
    medians_list = c()
    for (col in 1:ncol(normalization_factor_matrix)){
        medians_list = c(medians_list, median(normalization_factor_matrix[,col]))
        temp[,col] = temp[,col]/medians_list[col]
    }
    return(temp)
}

col_multiply = function(eset){
    result = matrix(data = 1,ncol=1,nrow=nrow(eset))
    for (column in 1:ncol(eset)){
        result = result*eset[:,column]
    }
    result = result^(1/ncol(eset))
    return(result)
}

normalization_methods$PoissonSeq = function(expression_set,subject_list){
    data_set = list(n=expression_set,y=subject_list)
    temp = expression_set
    normalization_factor_list = PoissonSeq::PS.Est.Depth(data_set$n)
    for (col in 1:ncol(expression_set)){
        temp[,col] = temp[,col]/normalization_factor_list[col]
    }
    return(temp)
}

# Simulating RNA-seq count matrix over 3 different datasets for 10k,20k and 50k with 8,60 and 550 subjects
simulation_options = list('sim1_opt' = RNAseq.SimOptions.2grp(ngenes=10000),
                          'sim2_opt' = RNAseq.SimOptions.2grp(ngenes=20000),
                          'sim3_opt' = RNAseq.SimOptions.2grp(ngenes=50000))

simulation_datasets = list('data1' = simRNAseq(simulation_options$sim1_opt,n1=4,n2=4),
                           'data2' = simRNAseq(simulation_options$sim2_opt,n1=4,n2=4),
                           'data3' = simRNAseq(simulation_options$sim3_opt,n1=4,n2=4),
                           'data4' = simRNAseq(simulation_options$sim1_opt,n1=50,n2=10),
                           'data5' = simRNAseq(simulation_options$sim2_opt,n1=50,n2=10),
                           'data6' = simRNAseq(simulation_options$sim3_opt,n1=50,n2=10),
                           'data7' = simRNAseq(simulation_options$sim1_opt,n1=450,n2=100),
                           'data8' = simRNAseq(simulation_options$sim2_opt,n1=450,n2=100),
                           'data9' = simRNAseq(simulation_options$sim3_opt,n1=450,n2=100))

# QC
# 1.Memory
# 2.Time
# 3.PCA

# 1.Memory graph and Time graph
# 2.Heatmap of genes - 2 main dataset
# 3.


