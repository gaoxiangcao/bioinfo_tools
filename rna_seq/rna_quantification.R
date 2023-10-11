####################################################################################################
# Data 2020-09-20
# Author Gaoxiang Cao
# E-mail gaoxiang_cao@163.com
####################################################################################################
rm(list = ls())

# make count table 
setwd("E:/RNA_Seq/data_process/DESeq2")
raw_df <- read.table(file = ".//out_table/293T-RNASeq-Ctrl_vs_KD.STAR.hg38.featureCounts.tsv",header = T,skip = 1,sep = "\t")
count_df <- raw_df[,c(7:10)]
rownames(count_df) <- as.character(raw_df[,1])
colnames(count_df) <- c("Ctrl_rep1","Ctrl_rep2","METTL3_KD_rep1","METTL3_KD_rep2")

# write.table(count_df, file = "./03.code_and_data/out_table/293T-RNASeq-Ctrl_vs_KD.STAR.hg38.featureCounts.FixColName.tsv",col.names = T,row.names = T,quote = F,sep = "\t")

###############################################################################
# Part I DESeq2
###############################################################################
# install DESeq2
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

# preparation 
rm(list = ls())
library(DESeq2)

# -------------------------------------------------------->>>>>>
# make obj
# -------------------------------------------------------->>>>>>

count_df <- read.table("./out_table/293T-RNASeq-Ctrl_vs_KD.STAR.hg38.featureCounts.FixColName.tsv", sep = "\t", header = T)
colnames(count_df)

# filter
count_df.filter <- count_df[rowSums(count_df) > 20 & apply(count_df, 1, function(x){all(x > 0) }),]
head(count_df.filter)

# condition table 
sample_df <- data.frame(
  condition = c(rep("ctrl",2), rep("KD",2)),
  cell_line = "293T"
)

# -------------------------------------------------------->>>>>>
# direct get result
# -------------------------------------------------------->>>>>>
rownames(sample_df) <- colnames(count_df.filter)

deseq2.obj <- DESeqDataSetFromMatrix(countData = count_df.filter, colData = sample_df, design = ~condition)
deseq2.obj

# test
deseq2.obj <- DESeq(deseq2.obj)

# get result
deseq2.obj.res <- results(deseq2.obj)
deseq2.obj.res.df <- as.data.frame(deseq2.obj.res)
head(deseq2.obj.res.df)

# -------------------------------------------------------->>>>>>>>>>
# step by step get test result 
# -------------------------------------------------------->>>>>>>>>>
eseq2.obj <- DESeqDataSetFromMatrix(countData = count_df, colData = sample_df, design = ~condition)
deseq2.obj

# normalization 
deseq2.obj <- estimateSizeFactors(deseq2.obj)
sizeFactors(deseq2.obj)

# dispersion
deseq2.obj <- estimateDispersions(deseq2.obj)
dispersions(deseq2.obj)

# plot dispersion
plotDispEsts(deseq2.obj, ymin = 1e-4)

# test 
deseq2.obj <- nbinomWaldTest(deseq2.obj)
deseq2.obj.res <- results(deseq2.obj)

# -------------------------------------------------------->>>>>>>>>>
# MA plot
# -------------------------------------------------------->>>>>>>>>>
DESeq2::plotMA(deseq2.obj.res, alpha=0.001)

# -------------------------------------------------------->>>>>>>>>>
# calculate with float
# -------------------------------------------------------->>>>>>>>>>
# count table 
rm(list=ls())
count_df <- read.table("./out_table/293T-RNASeq-Ctrl_vs_KD.STAR.hg38.featureCounts.FixColName.tsv", sep = "\t", header = T)

# filter 
colnames(count_df)
count_df.filter <- count_df[rowSums(count_df) > 20 & apply(count_df,1,function(x){ all(x > 0) }),]

# add decimal
count_df.filter.dec <- count_df.filter + abs(rnorm(nrow(count_df.filter) * ncol(count_df.filter)))

# condition table
sample_df <- data.frame(
  condition = c(rep("ctrl",2), rep("KD",2)),
  cell_line = "293T"
)
rownames(sample_df) <- colnames(count_df.filter.dec)

# ceiling floor round 
deseq2.obj.dec <- DESeqDataSetFromMatrix(countData = round(count_df.filter.dec,0), colData = sample_df, design = ~condition)
deseq2.obj.dec


a <- c(1:10)
b <- c(2:11)
c <- data.frame(a,b)
c <- c + abs(rnorm(nrow(c) * ncol(c)))
c

###############################################################################
# Part II edgeR
###############################################################################
rm(list=ls())

library(edgeR)

# -------------------------------------------------------->>>>>>>>>>
# make obj 
# -------------------------------------------------------->>>>>>>>>>
# count table 
count_df <- read.table(file = "./out_table/293T-RNASeq-Ctrl_vs_KD.STAR.hg38.featureCounts.FixColName.tsv",header = T,sep = "\t")

# filter 
colnames(count_df)
count_df.filter <- count_df[rowSums(count_df) > 20 & apply(count_df,1,function(x){ all(x > 0) }),]

# condition table
group_info = c(rep("ctrl",2), rep("KD",2))

dge.list.obj <- DGEList(counts = count_df.filter, group = group_info)
dge.list.obj

# -------------------------------------------------------->>>>>>>>>>
# Normalization
# -------------------------------------------------------->>>>>>>>>>
# Normalization method: "TMM","TMMwsp","RLE","upperquartile","none"
dge.list.obj <- calcNormFactors(dge.list.obj,method = "TMM")
dge.list.obj$samples

dge.list.obj <- calcNormFactors(dge.list.obj,method = "TMMwsp")
dge.list.obj$samples

dge.list.obj <- calcNormFactors(dge.list.obj,method = "upperquartile")
dge.list.obj$samples

dge.list.obj <- calcNormFactors(dge.list.obj,method = "RLE") # DESeq2, cuffdiff
dge.list.obj$samples

dge.list.obj <- calcNormFactors(dge.list.obj,method = "none")
dge.list.obj$samples

# raw data plot MDS
plotMDS(dge.list.obj)

# -------------------------------------------------------->>>>>>>>>>
# make design matrix
# -------------------------------------------------------->>>>>>>>>>
design.mat <- model.matrix(~group_info)

# -------------------------------------------------------->>>>>>>>>>
# estimate dispersion
# -------------------------------------------------------->>>>>>>>>>
dge.list.obj <- estimateDisp(dge.list.obj,design.mat)
dge.list.obj$common.dispersion
dge.list.obj$tagwise.dispersion

# 1st common dispersion
dge.list.obj <- estimateCommonDisp(dge.list.obj)

# 2nd tagwise dispersion
dge.list.obj <- estimateTagwiseDisp(dge.list.obj)

# plot dispersion
plotBCV(dge.list.obj, cex = 0.8)

# plot var and mean
plotMeanVar(dge.list.obj, show.raw=TRUE, show.tagwise=TRUE, show.binned=TRUE)

# -------------------------------------------------------->>>>>>>>>>
# test with likelihood ratio test
# -------------------------------------------------------->>>>>>>>>>
dge.list.res <- exactTest(dge.list.obj)
DEGs.res <- as.data.frame(topTags(dge.list.res,n=nrow(count_df.filter),sort.by = "logFC"))

# MA plot
select.sign.gene = decideTestsDGE(dge.list.res, p.value = 0.001) 
select.sign.gene_id = rownames(dge.list.res)[as.logical(select.sign.gene)]
plotSmear(dge.list.res, de.tags = select.sign.gene_id, cex = 0.5,ylim=c(-4,4)) 
abline(h = c(-2, 2), col = "blue")

# -------------------------------------------------------->>>>>>>>>>
# test with likelihood ratio test
# -------------------------------------------------------->>>>>>>>>>
fit <- glmFit(dge.list.obj, design.mat)
lrt <- glmLRT(fit, coef=2)
DEGs.res.lrt <- as.data.frame(topTags(lrt,n=nrow(count_df.filter),sort.by = "logFC"))

###############################################################################
# Part III limma
###############################################################################
rm(list=ls())

library(limma)

# -------------------------------------------------------->>>>>>>>>>
# make obj 
# -------------------------------------------------------->>>>>>>>>>
# count table 
count_df <- read.table(file = "./out_table/293T-RNASeq-Ctrl_vs_KD.STAR.hg38.featureCounts.FixColName.tsv",header = T,sep = "\t")

# filter 
colnames(count_df)
count_df.filter <- count_df[rowSums(count_df) > 20 & apply(count_df,1,function(x){ all(x > 0) }),]

# condition table
group_info = c(rep("ctrl",2), rep("KD",2))

design.mat <- model.matrix(~factor(group_info))

# normalization method
# “none”, “scale”, “quantile”, “Aquantile”, “Gquantile”, “Rquantile”, “Tquantile”, “vsn”, “cyclicloess”

voom.res <- limma::voom(count_df.filter, design = design.mat, normalize.method = "quantile")
voom.fit <- lmFit(voom.res, design.mat)
voom.eBayes <- eBayes(voom.fit)
voom.DEG <- topTable(voom.eBayes,coef = 2, number = Inf)

###############################################################################
# Part IV cuffdiff
###############################################################################
rm(list = ls())

setwd("E:/RNA_Seq/data_process/cuffdiff")

library(tidyverse)

# load table
cuffdiff_res <- read_tsv("./cuffdiff_result/gene_exp.diff")
colnames(cuffdiff_res)
colnames(cuffdiff_res)[10] <- "log2FC"

# filter 
cuffdiff_res.filter <- filter(cuffdiff_res, status == "OK")

# real sign
real_sign = rep("no",nrow(cuffdiff_res.filter))

select.FPKM <- (cuffdiff_res.filter$value_1 > 1 | cuffdiff_res.filter$value_2 > 1)
table(select.FPKM)

select.log2FC <- abs(cuffdiff_res.filter$log2FC) > 1
table(select.log2FC)

select.qval <- (cuffdiff_res.filter$q_value < 0.05)
table(select.qval)

real_sign[select.FPKM & select.log2FC & select.qval] <- "yes"
table(real_sign)

# MA plot
cuffdiff_res.filter <- mutate(cuffdiff_res.filter,
                              M.value = log2(value_2) - log2(value_1),
                              A.value = (log2(value_2) + log2(value_1)) / 2,
                              real_sign = real_sign)
cuffdiff_res.filter

# default
ggplot(cuffdiff_res.filter, aes(x=A.value, y=M.value)) + 
  geom_point(aes(color=significant)) + 
  geom_abline(intercept = 0,slope = 0,linetype="dashed") + 
  theme_bw() + 
  ylim(-4,4) + 
  scale_color_manual(values=c("gray", "red"))

# real sign 
ggplot(cuffdiff_res.filter, aes(x=A.value, y=M.value)) + 
  geom_point(aes(color=real_sign)) + 
  geom_abline(intercept = 0,slope = 0,linetype="dashed") + 
  theme_bw() + 
  ylim(-4,4) + 
  scale_color_manual(values=c("gray", "red"))

# volcano plot
ggplot(cuffdiff_res.filter, aes(x=log2FC, y=log10(q_value) * -1 )) + 
  geom_point(aes(color=significant)) + 
  geom_abline(intercept = 0,slope = 0,linetype="dashed") + 
  theme_bw() + 
  ylim(0,4) + 
  xlim(-5,5) + 
  scale_color_manual(values=c("gray", "red"))

