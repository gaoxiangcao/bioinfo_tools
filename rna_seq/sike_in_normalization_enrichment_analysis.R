# ---------------------------------------------------------------------->>>>>>>
# spike-in or housekeeping gene normalization
# ---------------------------------------------------------------------->>>>>>>
rm(list = ls())
library(RUVSeq)
setwd("E:/RNA_Seq/data_process")

# load raw count table 
count_df <- read.table(file = "./out_table/293T-RNASeq-Ctrl_vs_KD.STAR.hg38.featureCounts.FixColName.tsv", header = T, sep = "\t")
count_df

# load housekeeping gene list
hsk_gene_df <- read.table(file = "./data/human_housekeeping_gene_list.txt", header = T, sep = "\t")
hsk_gene_df

# check in list gene number
table(rownames(count_df) %in% hsk_gene_df$gene_id)
table(hsk_gene_df$gene_id %in% rownames(count_df))

# select Housekeeping gene list
hsk_gene_id.count <- count_df[rownames(count_df) %in% hsk_gene_df$gene_id, ]
hsk_gene_id.list <- rownames(count_df)[rownames(count_df) %in% hsk_gene_df$gene_id]

# select upper 100~200
mean.value <- apply(hsk_gene_id.count, 1, mean)
hsk_gene_id.list.select <- hsk_gene_id.list[order(mean.value, decreasing = T)][800:1000]

# make sample info table
sample_df <- data.frame(condition = c("ctrl", "ctrl", "KD", "KD"),
                        row.names = colnames(count_df))

# make obj
set.obj <- newSeqExpressionSet(counts = as.matrix(count_df), phenoData = sample_df)
plotPCA(set.obj)
plotPCA(set.obj, ylim = c(-1, 1))

# normalize by housekeeping gene
set.norm.obj <- RUVg(set.obj, hsk_gene_id.list.select, k = 2)
plotPCA(set.norm.obj)
plotPCA(set.norm.obj, ylim = c(-1, 1))

# check normal count
count_df.norm_hsk <- normCounts(set.norm.obj)
count_df.norm_hsk

# check normalize colSums count
colSums(count_df)
colSums(count_df.norm_hsk)

# ---------------------------------------------------->>>>>>>
# find DEGs with edgeR
# ---------------------------------------------------->>>>>>>
library(edgeR)

# condition table
group_info <- c(rep("ctrl", 2), rep("ko", 2))
group_info

# build obj
dge.list.obj <- DGEList(counts = normCounts(set.norm.obj), group = group_info)
dge.list.obj

# Norm -> estimate dispersion -> find DEGs
dge.list.obj <- calcNormFactors(dge.list.obj, method = "upperquartile")
dge.list.obj$samples

design.mat <- model.matrix(~group_info)
dge.list.obj <- estimateDisp(dge.list.obj, design.mat)
dge.list.obj$common.dispersion
dge.list.obj$tagwise.dispersion
# find DEGs
fit <- glmFit(dge.list.obj, design.mat)
lrt <- glmLRT(fit, coef=2)
DEGs.res.lrt <- as.data.frame(topTags(lrt,n=nrow(count_df.norm_hsk),sort.by = "logFC"))


# ---------------------------------------------------------------------->>>>>>>
# GO analysis
# ---------------------------------------------------------------------->>>>>>>
rm(list = ls())

library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)

# load table 
cuffdiff.res <- read.table("./cuffdiff_result/gene_exp.diff", header = T, sep = "\t")
cuffdiff.res <- dplyr::filter(cuffdiff.res, status == "OK")
colnames(cuffdiff.res)[10] <- "log2FC"
head(cuffdiff.res)

# filter
cuffdiff_res.filter <- filter(cuffdiff.res, status == "OK")

# 1.FPKM > 1
select.FPKM <- (cuffdiff.res$value_1 > 1 | cuffdiff.res$value_2 > 1)
table(select.FPKM)

# 2.log2FC > 1 or log2FC < 1
select.log2FC <- abs(cuffdiff.res$log2FC) > 1
table(select.log2FC)

# 3.FDR or p.adjust < 0.05
select.qval <- (cuffdiff.res$q_value < 0.05)
table(select.qval)

# select DEGs
select_vec <- (select.FPKM & select.log2FC & select.qval)
table(select_vec)

cuffdiff_res.filter.DEG <- as.vector(cuffdiff.res$gene_id[select_vec])
cuffdiff_res.filter.DEG

# GO  BP
enrich.go.BP <- enrichGO(
  gene = cuffdiff_res.filter.DEG, 
  OrgDb = org.Hs.eg.db, 
  keyType = "SYMBOL", 
  ont = "BP",
  pvalueCutoff = 0.5,
  qvalueCutoff = 0.5
  )

columns(org.Hs.eg.db)  # type of enrichment analysis
enrich.go.BP.df <- as.data.frame(enrich.go.BP)


barplot(enrich.go.BP)
dotplot(enrich.go.BP)

# ---------------------------------------------------------------------->>>>>>>
# KEGG analysis
# ---------------------------------------------------------------------->>>>>>>
# convert id 

DEG.entrez_id <- mapIds(x=org.Hs.eg.db, 
                        keys = cuffdiff_res.filter.DEG, 
                        keytype = "SYMBOL", 
                        column = "ENTREZID", 
                        multiVals = "first")
enrich.kegg.res <- enrichKEGG(DEG.entrez_id, organism = "hsa", keyType = "kegg")
barplot(enrich.kegg.res)
dotplot(enrich.kegg.res)

# syn test
syn_DEP_list_table <- readxl::read_excel("./test.xlsx", col_names = F)
syn_DEP_list_table <- syn_DEP_list_table %>% dplyr::filter(!str_detect(protein, "TRUE | ^S[A-Za-z0-9]*"))
colnames(syn_DEP_list_table) <- "protein"
syn_DEP_list <- as.vector(syn_DEP_list_table$protein)

# filter
enrich.kegg.res <- enrichKEGG(syn_DEP_list, 
                              organism = "syn", 
                              keyType = "kegg",
                              pvalueCutoff = 0.9
                              )
barplot(enrich.kegg.res)
dotplot(enrich.kegg.res)

# ---------------------------------------------------------------------->>>>>>>
# GSEA
# ---------------------------------------------------------------------->>>>>>>
# BiocManager::install("clusterProfiler")
# BiocManager::install("fgsea")
library(fgsea)
library(clusterProfiler)

# gene set and value infor
# run gsea

# load gmt
gene_list <- read.gmt(gmtfile = "./gesa/CELL_MIGRATION.geneset.gmt")
gene_list

# gene
fgsea.gene_list <- list()
fgsea.gene_list[["CELL_MIGRATION"]] <- as.character(gene_list$gene)

# value
fgsea.value <- cuffdiff.res$log2FC
names(fgsea.value) <- as.character(cuffdiff.res$gene_id)
fgsea.value.filter <- fgsea.value[!is.infinite(fgsea.value)]
fgsea.value.filter

# run 
fgsea.res <- fgsea(
  pathways = fgsea.gene_list, 
  stats = fgsea.value.filter,
  nperm = 1
  )
plotEnrichment(fgsea.gene_list[["CELL_MIGRATION"]], fgsea.value.filter)


# clusterProfiler gesa 
library(enrichplot)
data_sort <- cuffdiff.res %>% arrange(desc(log2FC))
data_sort <- data_sort %>% dplyr::filter((log2FC != Inf) & log2FC != -Inf)

gene_list <- data_sort$log2FC
names(gene_list) <- data_sort$gene_id
head(gene_list)

res <- gseGO(
  gene_list,    
  ont = "BP",    
  OrgDb = org.Hs.eg.db,   
  keyType = "SYMBOL",    
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",   
)

set_res <- res@result

gene_set <- read.gmt(gmtfile = "./gesa/CELL_MIGRATION.geneset.gmt")


gseaplot2(res, geneSetID = 1:3)
gseaplot2(res, geneSetID = 1:3)
# ---------------------------------------------------------------------->>>>>>>
# download org.Db
# ---------------------------------------------------------------------->>>>>>>
# install 
# BiocManager::install("AnnotationHub")
library(AnnotationHub)
hub <- AnnotationHub::AnnotationHub()

# query name 
query(hub, "Synechocystis")
query(hub, "Solanum")
query(hub, "Escherichia coli")

# download
Ecoli.OrgDb <- hub[["AH111572"]]
syn.OrgDb <- hub[["AH10608", force=T]]  # force download 
syn.OrgDb <- loadDb("./orgdb/hom.Synechocystis_sp..inp8.sqlite")
head(keys(syn.OrgDb))

syn_6803.OrgDb <- hub[["AH12836"]]

# download from https://annotationhub.bioconductor.org/
# load local file
library(AnnotationDbi)
syn_6803.db <- AnnotationDbi::loadDb("./orgdb/org.Synechocystis_sp._PCC_6803.eg.sqlite")
columns(syn_6803.db)
syn_6803.db

head(keys(syn_6803.db, keytype = "SYMBOL"))
head(keys(syn_6803.db, keytype = "ENTREZID"))
head(keys(syn_6803.db, keytype = "SYMBOL"))



