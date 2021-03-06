% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/4_multi_deg_all.R
\name{multi_deg_all}
\alias{multi_deg_all}
\title{multi_deg_all}
\usage{
multi_deg_all(
  exp,
  group_list,
  ids,
  logFC_cutoff = 1,
  pvalue_cutoff = 0.05,
  adjust = F,
  entriz = T,
  scale_before = F,
  n_cutoff = 3,
  cluster_cols = T,
  annotation_legend = F,
  show_rownames = F,
  legend = F,
  lab = NA,
  pkg = 4,
  symmetry = F,
  heat_union = T,
  heat_id = 1,
  gene_number = 200,
  color_volcano = c("blue", "grey", "red")
)
}
\arguments{
\item{exp}{A numeric matrix}

\item{group_list}{A factor with duplicated character or factor}

\item{ids}{a data.frame with 2 columns,including probe_id and symbol}

\item{logFC_cutoff}{Cutoff value of logFC,1 by defult.}

\item{pvalue_cutoff}{Cutoff value of pvalue,0.05 by defult.}

\item{adjust}{a logical value, would you like to use adjusted pvalue to draw this plot,FAlSE by defult.}

\item{scale_before}{deprecated parameter}

\item{n_cutoff}{3 by defalut , scale before plot and set a cutoff,usually 2 or 1.6}

\item{cluster_cols}{if F,heatmap will nor cluster in column}

\item{annotation_legend}{logical，show annotation legend or not}

\item{show_rownames}{logical,show rownames or not}

\item{legend}{logical,show legend or not}

\item{pkg}{a integer ,means which Differential analysis packages you used,we support three packages by now, 1,2,3,4 respectively means "DESeq2","edgeR","limma(voom)","limma"}

\item{symmetry}{a logical value ,would you like to get your plot symmetrical}

\item{heat_union}{logical ,use union or intersect DEGs for heatmap}

\item{heat_id}{id of heatmap,1 for all DEGs,2 for head and tail,3 for top n DEGs}

\item{gene_number}{how many DEGs will heatmap show .}
}
\value{
a list with deg data.frame, volcano plot and a list with DEGs.
}
\description{
do diffiencial analysis according to exprission set and group information
}
\examples{
gse = "GSE474"
geo = geo_download(gse)
geo$exp[1:4,1:4]
geo$exp=log2(geo$exp+1)
library(stringr)
group_list=ifelse(str_detect(geo$pd$title,"MObese"),"MObese",ifelse(str_detect(geo$pd$title,"NonObese"),"NonObese","Obese"))
group_list=factor(group_list,levels = c("NonObese","Obese","MObese"))
find_anno(geo$gpl,install = T)
ids <- toTable(hgu133aSYMBOL)
dcp = multi_deg_all(geo$exp,
group_list,ids,adjust = F)
dcp[[3]]
}
\seealso{
\code{\link{geo_download}};\code{\link{draw_volcano}};\code{\link{draw_venn}}
}
\author{
Xiaojie Sun
}
