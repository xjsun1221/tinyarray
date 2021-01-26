##' draw PCA plots
##'
##' do PCA analysis and print a PCA plot
##'
##' @param exp A numeric matrix
##' @param group_list A factor with duplicated character or factor
##' @param color color vector
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
##' group_list <- factor(rep(c("A","B"),each = 3))
##' draw_pca(exp,group_list)
##' draw_pca(exp,group_list,color = c("blue","red"))
##' @seealso
##' \code{\link{draw_heatmap}};\code{\link{draw_volcano}};\code{\link{draw_venn}}

draw_pca <-  function(exp,group_list,
                      color = c("#92C5DE","#F4A582","#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494","#B3B3B3")
){
  p1 <-  all(apply(exp,2,is.numeric))
  if(!p1) stop("exp must be a numeric matrix")
  p2  <-  (sum(!duplicated(group_list)) > 1)
  if(!p2) stop("group_list must more than 1")
  p3 <- is.factor(group_list)
  if(!p3) {
    group_list = factor(group_list)
    warning("group_list was covert to factor")
  }
  dat <- as.data.frame(t(exp))
  dat.pca <- PCA(dat, graph = FALSE)
  col = color[1:length(levels(group_list))]
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
##' print a heatmap plot for expression matrix and group by group_list paramter, exp will be scaled.
##'
##' @inheritParams draw_pca
##' @param scale_before deprecated parameter
##' @param n_cutoff 3 by defalut , scale before plot and set a cutoff,usually 2 or 1.6
##' @param annotation_legend logical，show annotation legend or not
##' @param cluster_cols if F,heatmap will nor cluster in column
##' @param color color for heatmap
##' @param color_an color for column annotation
##' @param legend logical,show legend or not
##' @param show_rownames logical,show rownames or not
##' @return a heatmap plot according to \code{exp} and grouped by \code{group}.
##' @author Xiaojie Sun
##' @importFrom pheatmap pheatmap
##' @importFrom ggplotify as.ggplot
##' @importFrom RColorBrewer brewer.pal
##' @export
##' @examples
##' #example data
##' exp = matrix(abs(rnorm(60,sd = 16)),nrow = 10)
##' exp[,1:3] <- exp[,1:3]+20
##' colnames(exp) <- paste0("sample",1:6)
##' rownames(exp) <- paste0("gene",1:10)
##' exp[1:4,1:4]
##' group_list = factor(rep(c("A","B"),each = 3))
##' draw_heatmap(exp,group_list)
##' #use iris
##' n = t(iris[,1:4]);colnames(n) = 1:150
##' group_list = iris$Species
##' draw_heatmap(n,group_list)
##' draw_heatmap(n,group_list)
##' draw_heatmap(n,group_list,color = colorRampPalette(c("green","black","red"))(100),
##'              color_an = c("red","blue","pink") )
##' @seealso
##' \code{\link{draw_pca}};\code{\link{draw_volcano}};\code{\link{draw_venn}}

draw_heatmap <-  function(n,
                          group_list,
                          scale_before = F,
                          n_cutoff = 3,
                          cluster_cols = T,
                          legend = F,
                          show_rownames = F,
                          annotation_legend=F,
                          color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
                          color_an = c("#92C5DE","#F4A582","#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494","#B3B3B3")){
  n = as.data.frame(n)
  if(scale_before) {
    message("scale_before parameter is deprecated")
    scale_before = F
  }
  p1 <-  all(apply(n,2,is.numeric))
  if(!p1) stop("n must be a numeric matrix")
  p2  <-  (sum(!duplicated(group_list)) > 1)
  if(!p2) stop("group_list must more than 1")
  p3 <- is.factor(group_list)
  if(!p3) {
    group_list = factor(group_list)
    warning("group_list was covert to factor")
  }
  annotation_col=data.frame(group=group_list)
  rownames(annotation_col)=colnames(n)
  col = color_an[1:length(levels(group_list))]
  ann_colors = list(group = col)
  names(ann_colors$group)=levels(group_list)

  if(!scale_before){
    p = as.ggplot(pheatmap(n,
                           show_colnames =F,
                           show_rownames = show_rownames,
                           scale = "row",
                           color = color,
                           annotation_col=annotation_col,
                           annotation_colors = ann_colors,
                           cluster_cols = cluster_cols,
                           breaks = seq(-n_cutoff,n_cutoff,length.out = length(color)),
                           legend = legend,
                           slient = T,
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
##' @param color color vector
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
##' @importFrom ggplot2 unit
##' @importFrom ggplot2 theme
##' @export
##' @examples
##' data("des")
##' head(deseq_data)
##' draw_volcano(deseq_data)
##' draw_volcano(deseq_data,pvalue_cutoff = 0.01,logFC_cutoff = 2)
##' draw_volcano(deseq_data,color = c("darkgreen", "darkgrey", "#B2182B"))

##' @seealso
##' \code{\link{draw_heatmap}};\code{\link{draw_pca}};\code{\link{draw_venn}}

draw_volcano <- function(deg,lab=NA,pvalue_cutoff = 0.05,logFC_cutoff= 1,pkg = 1,adjust = F,symmetry = F,color = c("blue", "grey","red")){
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
                      ' down, ',
                      nrow(dat[dat$change =='UP',]),
                      ' up'
  )
  p <- ggplot(data = dat,
              aes(x = logFC,
                  y = -log10(P.value))) +
    geom_point(alpha=0.4, size=1.75,
               aes(color=change)) +
    scale_color_manual(values=color)+
    geom_vline(xintercept=c(-logFC_cutoff,logFC_cutoff),lty=4,col="black",lwd=0.8) +
    geom_hline(yintercept = -log10(pvalue_cutoff),lty=4,col="black",lwd=0.8) +
    theme_bw()+
    labs(title=this_tile , x=lab, y="")+
    theme(plot.title = element_text(hjust = 0.5))
  ce = ceiling(max(abs(dat$logFC)))
  if(symmetry)p <- p+scale_x_continuous(limits = c(-ce,ce),expand = c(0,0))
  return(p)
}


##' draw a venn plot
##'
##' print a venn plot for deg result created by three packages
##'
##' @param x a list for plot
##' @param name main of the plot
##' @param color color vector
##' @return a venn plot according to \code{x}, \code{y} and.\code{z} named "name" paramter
##' @author Xiaojie Sun
##' @importFrom VennDiagram venn.diagram
##' @importFrom ggplotify as.ggplot
##' @importFrom cowplot as_grob
##' @export
##' @examples
##' x = list(Deseq2=sample(1:100,30),edgeR = sample(1:100,30),limma = sample(1:100,30))
##' draw_venn(x,"test")
##' draw_venn(x,"test",color = c("darkgreen", "darkblue", "#B2182B"))
##' @seealso
##' \code{\link{draw_pca}};\code{\link{draw_volcano}};\code{\link{draw_heatmap}}

draw_venn <- function(x,name,color = c("#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494","#B3B3B3")){
  if(as.numeric(dev.cur())!=1) graphics.off()
  if(!is.list(x)) stop("x must be a list")
  if(length(x)>7) stop("why do you give me so many elements to compare, I reject!")
  col = color[1:length(x)]
  p = venn.diagram(x = x,
                   imagetype ="png",
                   filename=NULL,
                   lwd=1,#圈线粗度
                   lty=1, #圈线类型
                   col=col,
                   fill=col,
                   cat.col=col,
                   cat.cex = 1,
                   cat.dist = -0.15,
                   rotation.degree = 0,
                   main = name,
                   main.cex = 1,
                   cex=1,
                   alpha = 0.1,
                   reverse=TRUE)
  p = as.ggplot(as_grob(p))
  file.remove(dir(pattern = ("^VennDiagram.*log$")))
  return(p)
}


##' draw boxplot for expression
##'
##' draw boxplot for expression
##'
##' @inheritParams draw_pca
##' @param method one of kruskal.test,aov,t.test and wilcox.test
##' @param width wdith of boxplot and error bar
##' @param sort whether the boxplot will be sorted
##' @param drop whether to discard insignificant values
##' @param pvalue_cutoff if drop = T，genes with p-values below the threshold will be drawn
##' @param xlab title of the x axis
##' @param ylab title of the y axis
##' @param grouplab title of group legend
##' @param p.label whether to show p vlaue in the plot
##' @param add_error_bar whether to add error bar
##' @return a boxplot according to \code{exp} and grouped by \code{group}.
##' @author Xiaojie Sun
##' @importFrom tidyr gather
##' @importFrom ggpubr stat_compare_means
##' @importFrom tibble rownames_to_column
##' @importFrom ggplot2 stat_boxplot
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 geom_boxplot
##' @importFrom ggplot2 theme
##' @importFrom ggplot2 aes
##' @importFrom ggplot2 labs
##' @importFrom ggplot2 scale_fill_manual
##' @importFrom ggplot2 element_text
##' @export
##' @examples
##' draw_boxplot(t(iris[,1:4]),iris$Species)
##' exp <-  matrix(rnorm(60),nrow = 10)
##' colnames(exp) <- paste0("sample",1:6)
##' rownames(exp) <- paste0("gene",1:10)
##' exp[,4:6] = exp[,4:6] +10
##' exp[1:4,1:4]
##' group_list <- factor(rep(c("A","B"),each = 3))
##' draw_boxplot(exp,group_list)
##' draw_boxplot(exp,group_list,color = c("grey","red"))
##' @seealso
##' \code{\link{draw_heatmap}};\code{\link{draw_volcano}};\code{\link{draw_venn}}

draw_boxplot = function(exp,group_list,
                        method = "kruskal.test",
                        sort = T,
                        drop = F,
                        width = 0.5,
                        pvalue_cutoff = 0.05,
                        xlab = "Gene",
                        ylab = "Expression",
                        grouplab = "Group",
                        p.label = F,
                        add_error_bar = F,
                        color = c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854",
                                 "#FFD92F", "#E5C494", "#B3B3B3")){
  p0 <-  all(apply(exp,2,is.numeric)) & (!is.null(rownames(exp)))
  if(!p0) stop("exp must be a numeric matrix with rownames")
  p1 = method %in% c("kruskal.test","aov","t.test","wilcox.test")
  if(!p1) stop("method should be one of kruskal.test,aov,t.test and wilcox.test")
  if(length(unique(group_list))>2 & method %in% c("t.test","wilcox.test")) stop("group_list parameter must have exactly 2 group")
  p2  <-  (sum(!duplicated(group_list)) > 1)
  if(!p2) stop("group_list must more than 1")
  p3 <- is.factor(group_list)
  if(!p3) {
    group_list = factor(group_list)
    warning("group_list was covert to factor")
  }
  if(method=="kruskal.test"){
    x = apply(exp, 1, function(x){
      kruskal.test(x~group_list)$p.value
    })
  }else if(method=="aov"){
    x = apply(exp, 1, function(x){
      summary(aov(x~group_list))[[1]]$`Pr(>F)`[1]
    })
  }else if(method=="t.test"){
    x = apply(exp, 1, function(x){
      t.test(x~group_list)$p.value
    })
  }else if(method=="wilcox.test"){
    x = apply(exp, 1, function(x){
      wilcox.test(x~group_list)$p.value
    })
  }
  if(sum(x< pvalue_cutoff)==0) {
    message("No rows below the threshold,plot all rows")
    drop = F
  }
  if(drop){
    x = x[ x < pvalue_cutoff]
    exp = exp[rownames(exp) %in% names(x),]
  }

  if(length(x)>40) message("seems to be too many rows")
  dat = rownames_to_column(as.data.frame(t(exp)),var = "sample")
  dat$group = group_list
  dat = gather(dat,
               rows,exp,-sample,-group)
  if(sort){
    dat$rows = factor(dat$rows,
                      levels = names(sort(x)),
                      ordered = T)
  }
  col = color[1:length(levels(group_list))]
  p = ggplot(dat,aes(rows,exp,fill = group))+
    geom_boxplot( width = width)+
    theme_bw()+
    theme(legend.position = "top")+
    labs(fill = grouplab,
         x = xlab,
         y = ylab)+
    scale_fill_manual(values = col)
  if(add_error_bar) p = stat_boxplot(geom ='errorbar', width = width)
  if(!p.label){
    p = p + stat_compare_means(aes(group = group,label = ..p.signif..),method = method)
  }else{
    p = p + stat_compare_means(aes(group = group,label = ..p.format..),method = method)
  }
  if(length(x)>10) p = p + theme(axis.text.x = element_text(angle=50,vjust = 0.5))
  return(p)
}
