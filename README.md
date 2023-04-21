# tinyarray

### Introduction

Hi, I'm Xiao Jie. This is an R package I wrote based on my own data analysis needs. I'm glad you found it. I will update some useful functions here on the public account "bioinfoplanet" and also do some other sharing.

###  Installation

#### 1.online

```
if(!require(tinyarray))install.packages("tinyarray")
if(!require(devtools))install.packages("devtools")
if(!require(tinyarray))devtools::install_github("xjsun1221/tinyarray",upgrade = FALSE,dependencies = TRUE)
```

#### 2.local

Click the green button "code" on this page, then click "Download ZIP" to download it to your working directory. Install it with `devtools::install_local("tinyarray-master.zip",upgrade = F,dependencies = T)`.

### functions

#### 1.basic

draw_heatmap(),draw_volcano(),draw_venn(),draw_boxplot(),draw_KM(),draw_venn(),risk_plot()

ggheat() is a function from the ggplot2 package that can be used to create heatmaps. It is still relatively immature and mainly used for aligning plots and collecting legends.

Something about ggheat():
https://mp.weixin.qq.com/s/WhsBf6QAhVXeXeScM59cSA

#### 2.Downstream Analysis of Gene Expression Array Data from GEO Database

geo_download(): Provide a GEO number and get back the expression matrix, clinical information table, and platform number used.

find_anno(): Look up the annotation of the array platform.

get_deg(): Provide the array expression matrix, grouping information, probe annotation and get back the differential analysis results.

multi_deg(): Differential analysis for multiple groups (up to 5).

If you want to do differential analysis and get the common figures in one step, you can use get_deg_all() and multi_deg_all(). This part mainly integrates and simplifies the differential analysis of GEOquery, Annoprobe, and limma.

quick_enrich(): Simple and intuitive enrichment analysis.

double_enrich(): Separate enrichment of up- and down-regulated genes, combined with plotting.

https://mp.weixin.qq.com/s/YQQoDsE5JaKpgFGlbEfQNg

https://mp.weixin.qq.com/s/j5IB_MQ0zeOCe1j_ahwtdQ

#### 3.Exploring Expression Matrices

make_tcga_group(): Quickly get the grouping according to the TCGA sample naming rules.

sam_filter(): Remove duplicate tumor samples in TCGA.

match_exp_cl(): Match TCGA expression matrix with clinical information.

trans_array(): Replace the row names of the matrix, such as replacing the probe names of the expression matrix with gene names.

trans_exp(): Convert TCGA or TCGA+GTeX data to gene IDs (old version, genecode v22 or v23)

trans_exp_new(): Convert TCGA or TCGA+GTeX data to gene IDs(new versions)

t_choose(): Do t-tests for individual genes in batches.

cor.full() and cor.one(): Calculate correlations between genes in batches.

#### 4.Survival Analysis and Visualization
 
point_cut(): Calculate the best cutoff point for survival analysis in batches.

surv_KM(): Do KM survival analysis in batches, supporting grouping with the best cutoff point.

surv_cox(): Do single factor Cox in batches, supporting grouping with the best cutoff point.

risk_plot(): Risk factor three-way linkage.

https://mp.weixin.qq.com/s/WYBhGxfGg6QFUPHFBashaA

exp_boxplot(): Draw T-N boxplot for the interested genes.

exp_surv(): Draw KM-plot for the interested genes.

box_surv(): Draw boxplot and KM-plot for the interested genes.

#### 5.Something about network graph

hypertest(): Do hypergeometric distribution test for mRNA and lncRNA in batches.

plcortest(): Do correlation test for mRNA and lncRNA in batches.

https://www.yuque.com/xiaojiewanglezenmofenshen/bsgk2d/dt0isp

interaction_to_edges(): Generate the connection table for the network graph based on the relationship table.

edges_to_nodes(): Generate the node table based on the connection table.

#### 6.Tricks

dumd(): Count how many values each column of the data frame has.

intersect_all(): Take the intersection of any number of vectors.

union_all(): Take the union of any number of vectors.

