##' draw PCA plots
##'
##' do PCA analysis and print a PCA plot
##'
##' @param exp A numeric matrix
##' @param group_list A factor with duplicated character or factor
##' @return a pca plot according to \code{exp} and grouped by \code{group}.
##' @author Xiaojie Sun
##' @importFrom FactoMineR PCA
##' @importFrom factoextra fviz_pca_ind
##' @export
##' @examples
##' draw_pca(t(iris[,1:4]),iris$Species)
##' exp <-  matrix(rnorm(60),nrow = 10)
##' colnames(exp) <- paste0("sample",1:6)
##' rownames(exp) <- paste0("gene",1:10)
##' exp[1:4,1:4]
##' group_list <- rep(c("A","B"),each = 3)
##' draw_pca(exp,group_list)
##' @seealso
##' \code{\link{draw_heatmap}};\code{\link{draw_volcano}};\code{\link{draw_venn}}

draw_pca <-  function(exp,group_list){
  p1 <-  all(apply(exp,2,is.numeric))
  if(!p1) stop("exp must be a numeric matrix")
  p2  <-  (sum(!duplicated(group_list)) > 1)
  if(!p2) stop("group_list must more than 1")
  p3 <- is.factor(group_list)
  if(!p3) stop("group_list must be a factor")
  dat <- as.data.frame(t(exp))
  dat.pca <- PCA(dat, graph = FALSE)
  colset = c("#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494","#B3B3B3")
  col = colset[1:length(levels(group_list))]
  fviz_pca_ind(dat.pca,
               geom.ind = "point",
               col.ind = group_list,
               palette = col,
               addEllipses = TRUE,
               legend.title = "Groups")
}

if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))



##' draw a heatmap plot
##'
##' print a heatmap plot for expression matrix and group by group_list paramter
##'
##' @inheritParams draw_pca
##' @param scale_before logical,if TRUE ,scale before plot,if FALSE,ues scale = "row" param form pheatmap
##' @param n_cutoff 2 by defalut , scale before plot and set a cutoff,usually 2 or 1.6
##' @param annotation_legend logical，show annotation legend or not
##' @param cluster_cols if F,heatmap will nor cluster in column
##' @return a heatmap plot according to \code{exp} and grouped by \code{group}.
##' @author Xiaojie Sun
##' @importFrom pheatmap pheatmap
##' @importFrom ggplotify as.ggplot
##' @export
##' @examples
##' #use your example data
##' exp = matrix(abs(rnorm(60,sd = 16)),nrow = 10)
##' exp[,1:3] <- exp[,1:3]+20
##' colnames(exp) <- paste0("sample",1:6)
##' rownames(exp) <- paste0("gene",1:10)
##' exp[1:4,1:4]
##' group_list = rep(c("A","B"),each = 3)
##' draw_heatmap(exp,group_list)
##' #use iris
##' n = t(iris[,1:4]);colnames(n) = 1:150
##' group_list = iris$Species
##' draw_heatmap(n,group_list)
##' @seealso
##' \code{\link{draw_pca}};\code{\link{draw_volcano}};\code{\link{draw_venn}}

draw_heatmap <-  function(n,group_list,scale_before = F,n_cutoff = 2,cluster_cols = T,annotation_legend=F){
  p1 <-  all(apply(n,2,is.numeric))
  if(!p1) stop("n must be a numeric matrix")
  p2  <-  (sum(!duplicated(group_list)) > 1)
  if(!p2) stop("group_list must more than 1")
  p3 <- is.factor(group_list)
  if(!p3) stop("group_list must be a factor")
  annotation_col=data.frame(group=group_list)
  rownames(annotation_col)=colnames(n)
  library(RColorBrewer)
  colset = c("#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494","#B3B3B3")
  col = colset[1:length(levels(group_list))]
  ann_colors = list(group = col)
  names(ann_colors$group)=levels(group_list)
  if(scale_before){
    n = t(scale(t(n)))
    n[n>n_cutoff] = n_cutoff
    n[n< -n_cutoff] = -n_cutoff
    p = as.ggplot(pheatmap(n,
                           show_colnames =F,
                           show_rownames = F,
                           annotation_col=annotation_col,
                           legend = F,
                           cluster_cols = cluster_cols,
                           annotation_legend = F,
                           annotation_colors = ann_colors,
                           annotation_names_col = F
    ))
  }
  if(!scale_before){
    p = as.ggplot(pheatmap(n,
                           show_colnames =F,
                           show_rownames = F,
                           scale = "row",
                           annotation_col=annotation_col,
                           annotation_colors = ann_colors,
                           cluster_cols = cluster_cols,
                           legend = F,
                           annotation_legend = annotation_legend,
                           annotation_names_col = F))
  }
  return(p)
}


##' draw a volcano plot
##'
##' print a volcano plot for Differential analysis result in data.frame fomat.
##'
##' @param deg a data.frame created by Differential analysis
##' @param pvalue_cutoff Cutoff value of pvalue,0.05 by defult.
##' @param logFC_cutoff Cutoff value of logFC,1 by defult.
##' @param pkg a integer ,means which Differential analysis packages you used,we support three packages by now, 1,2,3,4 respectively means "DESeq2","edgeR","limma(voom)","limma"
##' @param adjust a logical value, would you like to use adjusted pvalue to draw this plot,FAlSE by defult.
##' @param symmetry a logical value ,would you like to get your plot symmetrical
##' @return a volcano plot according to logFC and P.value(or adjust P.value)
##' @author Xiaojie Sun
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 aes
##' @importFrom ggplot2 geom_point
##' @importFrom ggplot2 scale_color_manual
##' @importFrom ggplot2 geom_vline
##' @importFrom ggplot2 geom_hline
##' @importFrom ggplot2 theme_bw
##' @importFrom ggplot2 labs
##' @importFrom ggplot2 theme
##' @importFrom ggplot2 element_text
##' @importFrom ggplot2 scale_x_continuous
##' @export
##' @examples
##' data("des")
##' head(deseq_data)
##' draw_volcano(deseq_data)
##' draw_volcano(deseq_data,pvalue_cutoff = 0.01,logFC_cutoff = 2)
##' @seealso
##' \code{\link{draw_heatmap}};\code{\link{draw_pca}};\code{\link{draw_venn}}

draw_volcano <- function(deg,lab=NA,pvalue_cutoff = 0.05,logFC_cutoff= 1,pkg = 1,adjust = F,symmetry = F){
  if(!is.data.frame(deg)) stop("deg must be a data.frame created by Differential analysis")
  if(pvalue_cutoff>0.1)warning("Your pvalue_cutoff seems too large")
  if(pvalue_cutoff>=1)stop("pvalue_cutoff will never larger than 1")
  if(!adjust){
    dat = switch(EXPR = pkg,
                 v1 = deg[,c(2,5)],
                 v2 = deg[,c(1,4)],
                 v3 = deg[,c(1,4)],
                 v4 = deg[,c(1,4)])
  }else{
    dat = switch(EXPR = pkg,
                 v1 = deg[,c(2,6)],
                 v2 = deg[,c(1,5)],
                 v3 = deg[,c(1,5)],
                 v4 = deg[,c(1,5)])
  }
  colnames(dat)[1:2]=c("logFC","P.value")
  #logFC_cutoff <- with(dat,mean(abs(logFC)) + 2*sd(abs(logFC)) )
  dat$change = with(dat,ifelse(logFC>logFC_cutoff & P.value<pvalue_cutoff , "UP",
                               ifelse(logFC< -logFC_cutoff & P.value<pvalue_cutoff , "DOWN","NOT")))
  if(is.na(lab)) lab = c("DESeq2","edgeR","limma(voom)","limma")[pkg]
  this_tile <- paste0(nrow(dat[dat$change =='DOWN',]),
                      'down gene,',
                      nrow(dat[dat$change =='UP',]),
                      'up gene'
  )
  p <- ggplot(data = dat,
              aes(x = logFC,
                  y = -log10(P.value))) +
    geom_point(alpha=0.4, size=3.5,
               aes(color=change)) +
    scale_color_manual(values=c("blue", "grey","red"))+
    geom_vline(xintercept=c(-logFC_cutoff,logFC_cutoff),lty=4,col="black",lwd=0.8) +
    geom_hline(yintercept = -log10(pvalue_cutoff),lty=4,col="black",lwd=0.8) +
    theme_bw()+
    labs(title=this_tile , x=lab, y="")+
    theme(plot.title = element_text(hjust = 0.5))
  ce = ceiling(max(abs(dat$logFC)))
  if(symmetry)p <- p+scale_x_continuous(limits = c(-ce,ce),expand = c(0,0))
  return(p)
}
