
# Importing packages required for all analysis
library(GEOquery)
library(NOISeq)
library(DESeq2)
library(PoissonSeq)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(ggrepel)
library(clusterProfiler)
library(DEGreport)
library(DOSE)
library(pathview)
library(tidyverse)
library(dplyr)
library(PROPER)
library(bench)

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
    normalization_factor_matrix[is.na(normalization_factor_matrix)] = 1
    normalization_factor_matrix[normalization_factor_matrix == Inf] = 0
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
        result = result*eset[,column]
    }
    result = result^(1/ncol(eset))
    return(result)
}

get_subject_list = function(expression_data,key){
    list = as.character(expression_data)
    subjects = grepl(key,list)
    return(as.integer(subjects))
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


# Running simulations on all three with benchmarking

analysis_DESeq2 = bench(

    DESeq2_data1 = normalization_methods$DESeq2(simulation_datasets$data1$counts),
    DESeq2_data2 = normalization_methods$DESeq2(simulation_datasets$data2$counts),
    DESeq2_data3 = normalization_methods$DESeq2(simulation_datasets$data3$counts),
    DESeq2_data4 = normalization_methods$DESeq2(simulation_datasets$data4$counts),
    DESeq2_data5 = normalization_methods$DESeq2(simulation_datasets$data5$counts),
    DESeq2_data6 = normalization_methods$DESeq2(simulation_datasets$data6$counts),
    DESeq2_data7 = normalization_methods$DESeq2(simulation_datasets$data7$counts),
    DESeq2_data8 = normalization_methods$DESeq2(simulation_datasets$data8$counts),
    DESeq2_data9 = normalization_methods$DESeq2(simulation_datasets$data9$counts),check = FALSE)
jpeg("DESeq2_Autoplot.jpeg")
autoplot(analysis_DESeq2)
dev.off()
write.table(analysis_DESeq2,"DESeq2_Memory_Time_Analysis.txt",quote=F,row.names=F)

analysis_PoissonSeq = bench(

    PoissonSeq_data1 = normalization_methods$PoissonSeq(simulation_datasets$data1$counts,simulation_datasets$data1$designs),
    PoissonSeq_data2 = normalization_methods$PoissonSeq(simulation_datasets$data2$counts,simulation_datasets$data2$designs),
    PoissonSeq_data3 = normalization_methods$PoissonSeq(simulation_datasets$data3$counts,simulation_datasets$data3$designs),
    PoissonSeq_data4 = normalization_methods$PoissonSeq(simulation_datasets$data4$counts,simulation_datasets$data4$designs),
    PoissonSeq_data5 = normalization_methods$PoissonSeq(simulation_datasets$data5$counts,simulation_datasets$data5$designs),
    PoissonSeq_data6 = normalization_methods$PoissonSeq(simulation_datasets$data6$counts,simulation_datasets$data6$designs),
    PoissonSeq_data7 = normalization_methods$PoissonSeq(simulation_datasets$data7$counts,simulation_datasets$data7$designs),
    PoissonSeq_data8 = normalization_methods$PoissonSeq(simulation_datasets$data8$counts,simulation_datasets$data8$designs),
    PoissonSeq_data9 = normalization_methods$PoissonSeq(simulation_datasets$data9$counts,simulation_datasets$data9$designs),check = FALSE)
jpeg("PoissonSeq_Autoplot.jpeg")
autoplot(analysis_PoissonSeq)
dev.off()
write.table(analysis_PoissonSeq,"PoissonSeq_Memory_Time_Analysis.txt",quote=F,row.names=F)

analysis_TMM = bench(

    TMM_data1 = normalization_methods$TMM(simulation_datasets$data1$counts),
    TMM_data2 = normalization_methods$TMM(simulation_datasets$data2$counts),
    TMM_data3 = normalization_methods$TMM(simulation_datasets$data3$counts),
    TMM_data4 = normalization_methods$TMM(simulation_datasets$data4$counts),
    TMM_data5 = normalization_methods$TMM(simulation_datasets$data5$counts),
    TMM_data6 = normalization_methods$TMM(simulation_datasets$data6$counts),
    TMM_data7 = normalization_methods$TMM(simulation_datasets$data7$counts),
    TMM_data8 = normalization_methods$TMM(simulation_datasets$data8$counts),
    TMM_data9 = normalization_methods$TMM(simulation_datasets$data9$counts),check = FALSE)
jpeg("TMM_Autoplot.jpeg")
autoplot(analysis_TMM)
dev.off()
write.table(analysis_TMM,"TMM_Memory_Time_Analysis.txt",quote=F,row.names=F)


analysis_RPKM = bench(

    RPKM_data1 = normalization_methods$RPKM(simulation_datasets$data1$counts),
    RPKM_data2 = normalization_methods$RPKM(simulation_datasets$data2$counts),
    RPKM_data3 = normalization_methods$RPKM(simulation_datasets$data3$counts),
    RPKM_data4 = normalization_methods$RPKM(simulation_datasets$data4$counts),
    RPKM_data5 = normalization_methods$RPKM(simulation_datasets$data5$counts),
    RPKM_data6 = normalization_methods$RPKM(simulation_datasets$data6$counts),
    RPKM_data7 = normalization_methods$RPKM(simulation_datasets$data7$counts),
    RPKM_data8 = normalization_methods$RPKM(simulation_datasets$data8$counts),
    RPKM_data9 = normalization_methods$RPKM(simulation_datasets$data9$counts),check = FALSE)
jpeg("RPKM_Autoplot.jpeg")
autoplot(analysis_RPKM)
dev.off()
write.table(analysis_RPKM,"RPKM_Memory_Time_Analysis.txt",quote=F,row.names=F)

analysis_UQ = bench(
    UQ_data1 = normalization_methods$UQ(simulation_datasets$data1$counts),
    UQ_data2 = normalization_methods$UQ(simulation_datasets$data2$counts),
    UQ_data3 = normalization_methods$UQ(simulation_datasets$data3$counts),
    UQ_data4 = normalization_methods$UQ(simulation_datasets$data4$counts),
    UQ_data5 = normalization_methods$UQ(simulation_datasets$data5$counts),
    UQ_data6 = normalization_methods$UQ(simulation_datasets$data6$counts),
    UQ_data7 = normalization_methods$UQ(simulation_datasets$data7$counts),
    UQ_data8 = normalization_methods$UQ(simulation_datasets$data8$counts),
    UQ_data9 = normalization_methods$UQ(simulation_datasets$data9$counts),check = FALSE)
jpeg("UQ_Autoplot.jpeg")
autoplot(analysis_UQ)
dev.off()
write.table(analysis_UQ,"UQ_Memory_Time_Analysis.txt",quote=F,row.names=F)

# Running the real data

analysis_GSE40562 = bench(
    DESeq2_R1 = normalization_methods$DESeq2(GSE40562_eset),
    PoissonSeq_R1 = normalization_methods$PoissonSeq(GSE40562_eset,get_subject_list(GSE40562$source_name_ch1,"normal")),
    TMM_R1 = normalization_methods$TMM(GSE40562_eset),
    RPKM_R1 = normalization_methods$RPKM(GSE40562_eset),
    UQ_R1 = normalization_methods$UQ(GSE40562_eset)
)
jpeg("GSE40562.jpeg")
autoplot(analysis_GSE40562)
dev.off()
write.table(analysis_GSE40562,"GSE40562_Memory_Time_Analysis.txt",quote=F,row.names=F)
DESeq2_R1 = normalization_methods$DESeq2(GSE40562_eset)
PoissonSeq_R1 = normalization_methods$PoissonSeq(GSE40562_eset,get_subject_list(GSE40562$source_name_ch1,"normal"))
TMM_R1 = normalization_methods$TMM(GSE40562_eset)
RPKM_R1 = normalization_methods$RPKM(GSE40562_eset)
UQ_R1 = normalization_methods$UQ(GSE40562_eset)
Raw_R1 = GSE40562_eset

analysis_GSE98582 = bench(
    DESeq2_R2 = normalization_methods$DESeq2(GSE98582_eset),
    PoissonSeq_R2 = normalization_methods$PoissonSeq(GSE98582_eset,get_subject_list(GSE98582$characteristics_ch1.2,"Control")),
    TMM_R2 = normalization_methods$TMM(GSE98582_eset),
    RPKM_R2 = normalization_methods$RPKM(GSE98582_eset),
    UQ_R2 = normalization_methods$UQ(GSE98582_eset)
)
jpeg("GSE98582.jpeg")
autoplot(analysis_GSE98582)
dev.off()
write.table(analysis_GSE98582,"GSE98582_Memory_Time_Analysis.txt",quote=F,row.names=F)
DESeq2_R2 = normalization_methods$DESeq2(GSE98582_eset)
PoissonSeq_R2 = normalization_methods$PoissonSeq(GSE98582_eset,get_subject_list(GSE98582$characteristics_ch1.2,"Control"))
TMM_R2 = normalization_methods$TMM(GSE98582_eset)
RPKM_R2 = normalization_methods$RPKM(GSE98582_eset)
UQ_R2 = normalization_methods$UQ(GSE98582_eset)
Raw_R2 = GSE98582_eset

# Comparing the post normalization count matrixes
# Boxplots
jpeg("GSE40562_BoxPlot.jpeg")
par(mfrow=c(2,3))
boxplot(log2(DESeq2_R1))
boxplot(log2(PoissonSeq_R1))
boxplot(log2(TMM_R1))
boxplot(log2(RPKM_R1))
boxplot(log2(UQ_R1))
boxplot(log2(Raw_R1))
dev.off()

jpeg("GSE98582_BoxPlot.jpeg")
par(mfrow=c(2,3))
boxplot(log2(DESeq2_R2))
boxplot(log2(PoissonSeq_R2))
boxplot(log2(TMM_R2))
boxplot(log2(RPKM_R2))
boxplot(log2(UQ_R2))
boxplot(log2(Raw_R2))
dev.off()

# Pheatmaps
jpeg("GSE40562_PheatMaps.jpeg")
par(mfrow=c(2,3))
pheatmap(log2(DESeq2_R1))
pheatmap(log2(PoissonSeq_R1))
pheatmap(log2(TMM_R1))
pheatmap(log2(RPKM_R1))
pheatmap(log2(UQ_R1))
pheatmap(log2(Raw_R1))
dev.off()

jpeg("GSE98582_PheatMap.jpeg")
par(mfrow=c(2,3))
pheatmap(log2(DESeq2_R2))
pheatmap(log2(PoissonSeq_R2))
pheatmap(log2(TMM_R2))
pheatmap(log2(RPKM_R2))
pheatmap(log2(UQ_R2))
pheatmap(log2(Raw_R2))
dev.off()


