# -------------------------------------------------------->>>>>>
# Part I cluster
# -------------------------------------------------------->>>>>>
rm(list = ls())
setwd("E:/RNA_Seq/data_process")

# -------------------------------------------------------->>>>>>>>>>
# load table and make matrix
# -------------------------------------------------------->>>>>>>>>>
single_cell.df <- read.table("./multiple_samples/GSE44183_human_expression_mat.txt", sep = "\t", header = T)
head(single_cell.df)
single_cell.mat <- as.matrix(single_cell.df[, -1])
rownames(single_cell.mat) <- as.character(single_cell.df$gene)

# -------------------------------------------------------->>>>>>>>>>
# run cluster
# -------------------------------------------------------->>>>>>>>>>
library(cluster)
dist.mat <- dist(t(single_cell.mat), method = "euclidean")
plot(hclust(dist.mat, method = "ward.D2"))

dist.mat <- dist(t(single_cell.mat), method = "manhattan")
plot(hclust(dist.mat))

hclust.res <- hclust(dist.mat)
hclust.res$height
hclust.res$labels

# Correlation 
cor.res <- cor(single_cell.mat, method = "spearman")
cor.res


library(pheatmap)
pheatmap(cor.res)

# set color
library(RColorBrewer)
RColorBrewer::display.brewer.all()
color.map = colorRampPalette(brewer.pal(n = 7, name = "Greens")) 
color.breaks = seq(0,1,length.out=101)
cor_matrix = cor(cor.res)
pheatmap(cor_matrix,color = color.map(100),breaks = color.breaks,show_rownames = T,main = "Pearson correlation and cluster plot")

# -------------------------------------------------------->>>>>>
# PCA
# -------------------------------------------------------->>>>>>

# PCA 
rm(list = ls())
setwd("E:/RNA_Seq/data_process")
single_cell.df <- read.table("./multiple_samples/GSE44183_human_expression_mat.txt", sep = "\t", header = T)
head(single_cell.df)
single_cell.mat <- as.matrix(single_cell.df[, -1])
rownames(single_cell.mat) <- as.character(single_cell.df$gene)

library(corrplot)
corrplot(cor(single_cell.mat))

# filter data
SD.vector = apply(single_cell.mat,1,FUN=function(x){sd(x)})
single_cell.mat.sort = single_cell.mat[order(SD.vector,decreasing = T),]
single_cell.mat.top3000 = single_cell.mat.sort[1:3000,]

# corrplot(cor(t(single_cell.mat.top3000)))

# PCA
library(corrplot)
library(FactoMineR)
library(factoextra)
library(RColorBrewer)

pca.res <- PCA(t(single_cell.mat.top3000), scale.unit = T, graph = F)

# show result
pca.res$var

# var info
summary(pca.res$var)

# plot scree
fviz_eig(pca.res, addlabels = T)

# plot effective factor
fviz_pca_var(pca.res, col.var = "coord")

#### plot directly
fviz_pca_ind(pca.res,label="all")

# -------------------------------------------------------->>>>>>
# MDA
# -------------------------------------------------------->>>>>>
# run MDS
dist.mat = dist(t(single_cell.mat), method = "minkowski", p=5)

MDS_res = cmdscale(dist.mat,eig=TRUE, k=2)

MDS_xy.df = data.frame(MDS_1 = MDS_res$points[,1],
                       MDS_2 = MDS_res$points[,2],
                       cell_info = colnames(single_cell.mat))

# plot with ggplot2
library(ggplot2)

ggplot(data = MDS_xy.df, aes(MDS_1,MDS_2,color=cell_info)) + 
  geom_point()  + 
  geom_text(aes(label=cell_info)) + 
  theme_bw()

# ---------------------------------------------------------------------------->>>>>>>>
# t-SNE
# ---------------------------------------------------------------------------->>>>>>>>

rm(list=ls())

setwd("E:/RNA_Seq/data_process")

# load package 
library(Rtsne)

# load table and make matrix
single_cell.df <- read.table(file = "./multiple_samples/GSE44183_human_expression_mat.txt",header = T,sep = "\t")

single_cell.mat <- as.matrix(single_cell.df[,-1])
rownames(single_cell.mat) <- as.character(single_cell.df$gene)

single_cell.mat.log2 = log2(single_cell.mat + 1)

SD.vector = apply(single_cell.mat.log2,1,FUN=function(x){sd(x)})

single_cell.mat.log2.sort = single_cell.mat.log2[order(SD.vector,decreasing = T),]
single_cell.mat.log2.top3000 = single_cell.mat.log2.sort[1:3000,]


set.seed(20200927)
tSNE_res = Rtsne(t(single_cell.mat.log2.top3000),dims = 3,perplexity = 2,pca = T)

# plot region 
tSNE_res.df = as.data.frame(tSNE_res$Y)
colnames(tSNE_res.df) = c("tSNE1","tSNE2","tSNE3")
tSNE_res.df$cell_info = colnames(single_cell.mat)

ggplot(data = tSNE_res.df, aes(tSNE1,tSNE2,color=cell_info)) + 
  geom_point()  + 
  geom_text(aes(label=cell_info)) + 
  theme_bw()

ggplot(data = tSNE_res.df, aes(tSNE1,tSNE3,color=cell_info)) + 
  geom_point()  + 
  geom_text(aes(label=cell_info)) + 
  theme_bw()

ggplot(data = tSNE_res.df, aes(tSNE2,tSNE3,color=cell_info)) + 
  geom_point()  + 
  geom_text(aes(label=cell_info)) + 
  theme_bw()

# use plotly
library(plotly)

p <- plot_ly(tSNE_res.df, x = ~tSNE1, y = ~tSNE2, z = ~tSNE3, color = ~cell_info) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'tSNE1'),
                      yaxis = list(title = 'tSNE2'),
                      zaxis = list(title = 'tSNE3')))
p















