##' draw PCA plots
##'
##' do PCA analysis and print a PCA plot
##'
##' @param exp A numeric matrix
##' @param group_list A factor with duplicated character or factor
##' @param color color vector
##' @param addEllipses logical,add ellipses or not
##' @param style plot style,"default","ggplot2"and "3D"
##' @param color.label color legend label
##' @param title plot title
##' @param ... other paramters from fviz_pca_ind
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 scale_color_manual
##' @importFrom ggplot2 scale_fill_manual
##' @importFrom ggplot2 theme_classic
##' @importFrom ggplot2 theme
##' @importFrom ggplot2 labs
##' @importFrom ggplot2 aes
##' @importFrom ggplot2 stat_ellipse
##' @return a pca plot according to \code{exp} and grouped by \code{group}.
##' @author Xiaojie Sun
##' @export
##' @examples
##' draw_pca(t(iris[,1:4]),iris$Species)
##' draw_pca(t(iris[,1:4]),iris$Species,style = "ggplot2")
##' draw_pca(t(iris[,1:4]),iris$Species,style = "3D")
##' #change color
##' draw_pca(t(iris[,1:4]),iris$Species,color = c("#E78AC3", "#A6D854", "#FFD92F"))


##' @seealso
##' \code{\link{draw_heatmap}};\code{\link{draw_volcano}};\code{\link{draw_venn}}

draw_pca <-  function(exp,group_list,
                      color = c("#2874C5","#f87669","#e6b707","#868686","#92C5DE","#F4A582","#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494","#B3B3B3"),
                      addEllipses = TRUE,
                      style = "default",
                      color.label = "Group",
                      title = "",
                      ...){
  if(!requireNamespace("FactoMineR",quietly = TRUE)) {
    stop("Package \"FactoMineR\" needed for this function to work.
         Please install it by install.packages('FactoMineR')",call. = FALSE)
  }
  if(!requireNamespace("factoextra",quietly = TRUE)) {
    stop("Package \"factoextra\" needed for this function to work.
         Please install it by install.packages('factoextra')",call. = FALSE)
  }
  p1 <-  all(apply(exp,2,is.numeric))
  if(!p1) stop("exp must be a numeric matrix")
  p2  <-  (sum(!duplicated(group_list)) > 1)
  if(!p2) stop("group_list must more than 1")
  dat <- as.data.frame(t(exp))
  dat.pca <- FactoMineR::PCA(dat, graph = FALSE)
  col = color[1:length(levels(group_list))]
  if(style == "default"){
    factoextra::fviz_pca_ind(dat.pca,
                             geom.ind = "point",
                             col.ind = group_list,
                             addEllipses = addEllipses,
                             palette = col,
                             legend.title = "Groups",
                             title = title,
                             ...)
  }else if(style == "ggplot2"){
    pdat = data.frame(dat.pca[["ind"]][["coord"]],
                      Group = group_list)
    p = ggplot(pdat,aes(Dim.1,Dim.2))+
      geom_point(aes(Dim.1,Dim.2,fill = Group),
                 shape = 21,color = "black")+
      scale_color_manual(values = color[1:nlevels(group_list)])+
      scale_fill_manual(values = color[1:nlevels(group_list)])+
      theme_classic()+
      theme(legend.position = "top")+
      labs(color = color.label,fill = color.label,title = title)

    if(addEllipses) p = p +
      stat_ellipse(aes(color = Group,fill = Group),
                   geom = "polygon",
                   alpha = 0.3,
                   linetype = 2)
    return(p)
  }else if(style == "3D"){
    colors = color[as.numeric(group_list)]
    pdat = data.frame(dat.pca[["ind"]][["coord"]],
                      Group = group_list)
    scatterplot3d::scatterplot3d(pdat[,1:3],
                                 color = "black",
                                 pch = 21,
                                 bg = colors,
                                 main = title)
    graphics::legend("bottom",col = "black",
           legend = levels(group_list),
           pt.bg =  color[1:nlevels(group_list)], pch = 21,
           inset = -0.2, xpd = TRUE,
           horiz = TRUE)
  }
}

if(getRversion() >= "2.15.1")  utils::globalVariables(c(".","Dim.1","Dim.2","Group","pdat"))



##' draw a heatmap plot
##'
##' print a heatmap plot for expression matrix and group by group_list praramter, exp will be scaled.
##'
##' @inheritParams draw_pca
##' @param n  A numeric matrix
##' @param scale_before deprecated parameter
##' @param n_cutoff 3 by defalut , scale before plot and set a cutoff,usually 2 or 1.6
##' @param annotation_legend logical,show annotation legend or not
##' @param color color for heatmap
##' @param color_an color for column annotation
##' @param legend logical,show legend or not
##' @param show_rownames logical,show rownames or not
##' @param scale logical,scale the matrix or not
##' @param main the title of the plot
##' @param split_column split column by group_list
##' @param show_column_title show column title or not
##' @param ... other parameters from pheatmap
##' @return a heatmap plot according to \code{exp} and grouped by \code{group}.
##' @author Xiaojie Sun
##' @importFrom pheatmap pheatmap
##' @export
##' @examples
##' #example data
##' exp = matrix(abs(rnorm(60,sd = 16)),nrow = 10)
##' exp[,4:6] <- exp[,4:6]+20
##' colnames(exp) <- paste0("sample",1:6)
##' rownames(exp) <- paste0("gene",1:10)
##' exp[1:4,1:4]
##' group_list = factor(rep(c("A","B"),each = 3))
##' draw_heatmap(exp,group_list)
##' #use iris
##' n = t(iris[,1:4]);colnames(n) = 1:150
##' group_list = iris$Species
##' draw_heatmap(n,group_list)
##' draw_heatmap(n,group_list,color = colorRampPalette(c("green","black","red"))(100),
##'              color_an = c("red","blue","pink") )
##' @seealso
##' \code{\link{draw_pca}};\code{\link{draw_volcano}};\code{\link{draw_venn}}

draw_heatmap <-  function(n,
                          group_list,
                          scale_before = FALSE,
                          n_cutoff = 3,
                          legend = FALSE,
                          show_rownames = FALSE,
                          annotation_legend=FALSE,
                          split_column = FALSE,
                          show_column_title = FALSE,
                          color = grDevices::colorRampPalette(c("#2fa1dd", "white", "#f87669"))(100),
                          color_an = c("#2fa1dd","#f87669","#e6b707","#868686","#92C5DE","#F4A582","#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494","#B3B3B3"),
                          scale = TRUE,
                          main = NA,...){
  if(!requireNamespace("ggplotify",quietly = TRUE)) {
    stop("Package \"ggplotify\" needed for this function to work. Please install it by install.packages('ggplotify')",call. = FALSE)
  }
  p3 <- is.factor(group_list)
  if(!p3) {
    group_list = factor(group_list)
    warning("group_list was covert to factor")
  }
  n = as.data.frame(n)
  if(scale_before) {
    message("scale_before parameter is deprecated")
    scale_before = FALSE
  }
  p1 <-  all(apply(n,2,is.numeric))
  if(!p1) stop("non-numeric matrix detected")
  p2  <-  (sum(!duplicated(group_list)) > 1)
  if(!p2) stop("group_list must more than 1")
  if(split_column){
    if(!requireNamespace("circlize",quietly = TRUE)) {
      stop("Package \"circlize\" needed for this function to work. Please install it by install.packages('circlize')",call. = FALSE)
    }
    if(!requireNamespace("ComplexHeatmap",quietly = TRUE)) {
      stop("Package \"ComplexHeatmap\" needed for this function to work. Please install it by BiocManager::install('ComplexHeatmap')",call. = FALSE)
    }
    col_fun = circlize::colorRamp2(c(-n_cutoff, 0, n_cutoff), c(color[1], "white", color[100]))
    if(scale){
      n = t(scale(t(n)))
      n[n > n_cutoff] = n_cutoff
      n[n < - n_cutoff] = -n_cutoff
    }
    m0 = ComplexHeatmap::Heatmap(n,
                    column_split = group_list)
    column_title = suppressWarnings(names(ComplexHeatmap::column_order(m0)))
    top_annotation = ComplexHeatmap::HeatmapAnnotation(
      cluster = ComplexHeatmap::anno_block(gp = grid::gpar(fill = (color_an[1:nlevels(group_list)])[match(column_title,levels(group_list))]),
                           labels = column_title,
                           labels_gp = grid::gpar(col = "white", fontsize = 12)))
    if(show_column_title){
      p = ComplexHeatmap::Heatmap(n,name = " ",
                                  col = col_fun,
                                  top_annotation = top_annotation,
                                  column_split = group_list,
                                  show_heatmap_legend = legend ,
                                  border = sum(dim(n))< 150,
                                  rect_gp = grid::gpar(col = "grey", lwd = 1),
                                  show_column_names = FALSE,
                                  show_row_names = show_rownames)
    }else{
      p = ComplexHeatmap::Heatmap(n,name = " ",
                                  col = col_fun,
                                  top_annotation = top_annotation,
                                  column_split = group_list,
                                  show_heatmap_legend = legend ,
                                  border = sum(dim(n))< 150,
                                  rect_gp = grid::gpar(col = "grey", lwd = 1),
                                  show_column_names = FALSE,
                                  show_row_names = show_rownames ,
                                  column_title = NULL)
    }
  }else{
    annotation_col=data.frame(group=group_list)
    rownames(annotation_col)=colnames(n)
    col = color_an[1:length(levels(group_list))]
    ann_colors = list(group = col)
    if(is.null(names(ann_colors$group)))names(ann_colors$group)=levels(group_list)
    if(scale) {
      scale_row = "row"
      breaks = seq(-n_cutoff,n_cutoff,length.out = length(color))
    } else {
      scale_row = "none"
      breaks = NA
    }
    if(!scale_before){
      p = pheatmap(n,
                   show_colnames =F,
                   show_rownames = show_rownames,
                   scale = scale_row,
                   color = color,
                   annotation_col=annotation_col,
                   annotation_colors = ann_colors,
                   breaks = breaks,
                   legend = legend,
                   silent = TRUE,
                   annotation_legend = annotation_legend,
                   main = main,...)
      p = ggplotify::as.ggplot(p)
    }
  }
  return(p)
}


##' draw a volcano plot
##'
##' print a volcano plot for Differential analysis result in data.frame format.
##'
##' @param deg a data.frame created by Differential analysis
##' @param pvalue_cutoff Cutoff value of pvalue,0.05 by default.
##' @param logFC_cutoff Cutoff value of logFC,1 by default.
##' @param pkg a integer ,means which Differential analysis packages you used,we support three packages by now, 1,2,3,4 respectively means "DESeq2","edgeR","limma(voom)","limma"
##' @param adjust a logical value, would you like to use adjusted pvalue to draw this plot,FAlSE by default.
##' @param symmetry a logical value ,would you like to get your plot symmetrical
##' @param color color vector
##' @param lab label for  x axis in volcano plot, if NA , x axis names by package
##' @param xlab.package whether to use the package name as the x axis name
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
##' head(deseq_data)
##' draw_volcano(deseq_data)
##' draw_volcano(deseq_data,pvalue_cutoff = 0.01,logFC_cutoff = 2)
##' draw_volcano(deseq_data,color = c("darkgreen", "darkgrey", "#B2182B"))

##' @seealso
##' \code{\link{draw_heatmap}};\code{\link{draw_pca}};\code{\link{draw_venn}}

draw_volcano <- function(deg,lab=NA,xlab.package = TRUE,
                         pvalue_cutoff = 0.05,logFC_cutoff= 1,
                         pkg = 1,adjust = FALSE,symmetry = FALSE,
                         color = c("#2874C5", "grey","#f87669")){
  if(!is.data.frame(deg)) stop("deg must be a data.frame created by Differential analysis")
  if(pvalue_cutoff>0.1)warning("Your pvalue_cutoff seems too large")
  if(pvalue_cutoff>=1)stop("pvalue_cutoff will never larger than 1")

  if(!adjust){
    dat = switch(EXPR = pkg,
                 v1 = deg[,c(2,5)],
                 v2 = deg[,c(1,4)],
                 v3 = deg[,c("logFC", "P.Value")],
                 v4 = deg[,c("logFC", "P.Value")])
  }else{
    dat = switch(EXPR = pkg,
                 v1 = deg[,c(2,6)],
                 v2 = deg[,c(1,5)],
                 v3 = deg[,c("logFC", "adj.P.Val")],
                 v4 = deg[,c("logFC", "adj.P.Val")])
  }
  colnames(dat)[1:2]=c("logFC","P.value")
  #logFC_cutoff <- with(dat,mean(abs(logFC)) + 2*sd(abs(logFC)) )
  if(is.null(deg$change)){

    dat$change = with(dat,ifelse(logFC>logFC_cutoff & P.value<pvalue_cutoff , "UP",
                                 ifelse(logFC< -logFC_cutoff & P.value<pvalue_cutoff , "DOWN","NOT")))
  }else{
    dat$change = str_to_upper(deg$change)
  }
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
    geom_vline(xintercept = c(-logFC_cutoff,logFC_cutoff),linetype=4,col="black",linewidth=0.8) +
    geom_hline(yintercept = -log10(pvalue_cutoff),linetype=4,col="black",linewidth=0.8) +
    theme_bw()+
    labs(title = this_tile , x = lab, y = "-log10(P.value)")+
    theme(plot.title = element_text(hjust = 0.5))
  ce = ceiling(max(abs(dat$logFC)))
  if(symmetry)p <- p+scale_x_continuous(limits = c(-ce,ce),expand = c(0,0))
  return(p)
}
utils::globalVariables(c("logFC","P.value","change"))

##' draw a venn plot
##'
##' print a venn plot for deg result created by three packages
##'
##' @inheritParams VennDiagram::venn.diagram
##' @inheritParams VennDiagram::draw.single.venn
##' @param x a list for plot
##' @param color color vector
##' @param reverse logical,reflect the three-set Venn diagram along its central
##' vertical axis of symmetry. Use in combination with rotation
##' to generate all possible set orders
##' @param ... other parameters from venn.diagram
##' @return a venn plot according to \code{x}, \code{y} and.\code{z} named "name" paramter
##' @author Xiaojie Sun
##' @export
##' @examples
##' x = list(Deseq2=sample(1:100,30),edgeR = sample(1:100,30),limma = sample(1:100,30))
##' draw_venn(x,"test")
##' draw_venn(x,"test",color = c("darkgreen", "darkblue", "#B2182B"))
##' @seealso
##' \code{\link{draw_pca}};\code{\link{draw_volcano}};\code{\link{draw_heatmap}}

draw_venn <- function(x,main,
                      color = c("#2874C5","#f87669","#e6b707","#868686","#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494","#B3B3B3"),
                      imagetype ="png",
                      filename=NULL,
                      lwd=1,
                      lty=1,
                      col=color[1:length(x)],
                      fill=color[1:length(x)],
                      cat.col=color[1:length(x)],
                      cat.cex = 1,
                      cat.dist = -0.15,
                      rotation.degree = 0,
                      main.cex = 1,
                      cex=1,
                      alpha = 0.1,
                      reverse=TRUE,
                      ...){
  if(as.numeric(grDevices::dev.cur())!=1) grDevices::graphics.off()
  if(!requireNamespace("VennDiagram",quietly = TRUE)) {
    stop("Package \"VennDiagram\" needed for this function to work. Please install it byby install.packages('VennDiagram')",call. = FALSE)
  }
  if(!requireNamespace("ggplotify",quietly = TRUE)) {
    stop("Package \"ggplotify\" needed for this function to work. Please install it by install.packages('ggplotify')",call. = FALSE)
  }
  if(!requireNamespace("cowplot",quietly = TRUE)) {
    stop("Package \"cowplot\" needed for this function to work. Please install it by install.packages('cowplot')",call. = FALSE)
  }
  if(!is.list(x)) stop("x must be a list")
  if(length(x)>7) stop("why do you give me so many elements to compare, I reject!")
  p = VennDiagram::venn.diagram(x = x,
                   imagetype =imagetype,
                   filename=filename,
                   lwd=lwd,
                   lty=lty,
                   col=col,
                   fill=fill,
                   cat.col=fill,
                   cat.cex = cat.cex,
                   cat.dist = cat.dist,
                   rotation.degree = rotation.degree,
                   main = main,
                   main.cex = main.cex,
                   cex=cex,
                   alpha = alpha,
                   reverse=reverse,
                   ...)
  p = ggplotify::as.ggplot(cowplot::as_grob(p))
  file.remove(dir(pattern = ("^VennDiagram.*log$")))
  return(p)
}


##' draw boxplot for expression
##'
##' draw boxplot for expression
##'
##' @inheritParams draw_pca
##' @param method one of kruskal.test,aov,t.test and wilcox.test
##' @param width width of boxplot and error bar
##' @param sort whether the boxplot will be sorted
##' @param drop whether to discard insignificant values
##' @param pvalue_cutoff if drop = TRUE,genes with p-values below the threshold will be drawn
##' @param xlab title of the x axis
##' @param ylab title of the y axis
##' @param grouplab title of group legend
##' @param p.label whether to show p value in the plot
##' @param add_error_bar whether to add error bar
##' @param ... other parameters from stat_compare_means
##' @return a boxplot according to \code{exp} and grouped by \code{group}.
##' @author Xiaojie Sun
##' @importFrom tibble rownames_to_column
##' @importFrom ggplot2 stat_boxplot
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 geom_boxplot
##' @importFrom ggplot2 theme
##' @importFrom ggplot2 aes
##' @importFrom ggplot2 labs
##' @importFrom ggplot2 scale_fill_manual
##' @importFrom ggplot2 element_text
##' @importFrom ggplot2 after_stat
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
                        sort = TRUE,
                        drop = FALSE,
                        width = 0.5,
                        pvalue_cutoff = 0.05,
                        xlab = "Gene",
                        ylab = "Expression",
                        grouplab = "Group",
                        p.label = FALSE,
                        add_error_bar = FALSE,
                        color = c("#2874C5","#f87669","#e6b707","#868686","#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854",
                                 "#FFD92F", "#E5C494", "#B3B3B3"),
                                 ...){
  if(!requireNamespace("ggpubr",quietly = TRUE)) {
    stop("Package \"ggpubr\" needed for this function to work.
         Please install it by install.packages('ggpubr')",call. = FALSE)
  }
  if(!requireNamespace("tidyr",quietly = TRUE)) {
    stop("Package \"tidyr\" needed for this function to work.
         Please install it by install.packages('tidyr')",call. = FALSE)
  }
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
      stats::kruskal.test(x~group_list)$p.value
    })
  }else if(method=="aov"){
    x = apply(exp, 1, function(x){
      summary(stats::aov(x~group_list))[[1]]$`Pr(>F)`[1]
    })
  }else if(method=="t.test"){
    x = apply(exp, 1, function(x){
      stats::t.test(x~group_list)$p.value
    })
  }else if(method=="wilcox.test"){
    x = apply(exp, 1, function(x){
      stats::wilcox.test(x~group_list)$p.value
    })
  }
  if(sum(x< pvalue_cutoff)==0) {
    message("No rows below the threshold,plot all rows")
    drop = FALSE
  }
  if(drop){
    x = x[ x < pvalue_cutoff]
    exp = exp[rownames(exp) %in% names(x),]
  }

  if(length(x)>40) message("it seems to be too many rows")
  dat = rownames_to_column(as.data.frame(t(exp)),var = "sample")
  dat$group = group_list
  dat = tidyr::gather(dat,
               rows,exp,-sample,-group)
  if(sort){
    dat$rows = factor(dat$rows,
                      levels = names(sort(x)),
                      ordered = TRUE)
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
    p = p + ggpubr::stat_compare_means(aes(group = group,
                                           label = ggplot2::after_stat(p.signif)),
                                       method = method,...)
  }else{
    p = p + ggpubr::stat_compare_means(aes(group = group,
                                           label = ggplot2::after_stat(p.format)),
                                       method = method,...)
  }
  if(length(x)>10) p = p +
    theme(axis.text.x = element_text(vjust = 1,hjust = 1,angle = 80))
  return(p)
}

utils::globalVariables(c(".","rows","group","p.signif","p.format"))

##' ggheat
##'
##' draw heatmap plot with annotation by ggplot2
##'
##' @param dat expression matrix for plot
##' @param group group for expression colnames
##' @param cluster logical,cluster in both rows and column or not, default F,now replaced by cluster_rows and cluster_cols.
##' @param cluster_rows logical, if rows (on the plot) should be clustered, default F
##' @param cluster_cols logical, if column (on the plot) should be clustered, default F
##' @param show_rownames logical,show rownames in plot or not, default T
##' @param show_colnames logical,show colnames in plot or not, default T
##' @param groupname name of group legend
##' @param expname name of exp legend
##' @param fill_mid use median value as geom_tile fill midpoint
##' @param color color for heatmap
##' @param legend_color color for legend
##' @return a ggplot object
##' @author Xiaojie Sun
##' @importFrom pheatmap pheatmap
##' @importFrom tibble rownames_to_column
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 geom_tile
##' @importFrom ggplot2 scale_fill_manual
##' @importFrom ggplot2 theme
##' @importFrom ggplot2 aes
##' @importFrom ggplot2 element_blank
##' @importFrom ggplot2 scale_x_continuous
##' @importFrom ggplot2 element_rect
##' @importFrom ggplot2 scale_fill_gradient2
##' @importFrom ggplot2 scale_x_discrete
##' @importFrom ggplot2 labs
##' @importFrom patchwork wrap_plots
##' @export
##' @examples
##' exp_dat = matrix(sample(100:1000,40),ncol = 4)
##' exp_dat[1:(nrow(exp_dat)/2),] =  exp_dat[1:(nrow(exp_dat)/2),]-1000
##' rownames(exp_dat) = paste0("sample",1:nrow(exp_dat))
##' colnames(exp_dat) = paste0("gene",1:ncol(exp_dat))

##' group = rep(c("A","B"),each = nrow(exp_dat)/2)
##' group = factor(group,levels = c("A","B"))
##' ggheat(exp_dat,group)
##' ggheat(exp_dat,group,cluster_rows = TRUE)
##' ggheat(exp_dat,group,cluster_rows = TRUE,show_rownames = FALSE,
##'        show_colnames = FALSE,groupname = "risk",expname = "expression")



ggheat = function(dat,group,cluster = FALSE,
                  color = c("#2874C5", "white","#f87669"),
                  legend_color = c("#2874C5","#f87669","#e6b707","#868686","#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F",
                                   "#E5C494", "#B3B3B3"),
                  show_rownames = TRUE,show_colnames = TRUE,
                  cluster_rows = FALSE,cluster_cols = FALSE,
                  groupname = "group",expname = "exp",
                  fill_mid = TRUE){
  dat = data.frame(dat)
  ph = pheatmap::pheatmap(t(dat),silent = TRUE)
  if(cluster){
    cluster_rows = TRUE
    cluster_cols = TRUE
  }
  if(cluster_rows){
    dat = dat[,ph$tree_row$order]
  }
  if(cluster_cols){
    dat = dat[ph$tree_col$order,]
    group = group[ph$tree_col$order]
  }

  dat$group = group
  dat2 = tidyr::gather(rownames_to_column(dat,
                                   "samples"),
                gene,exp,-group,-samples)

  dat2$samples = factor(dat2$samples,levels = rownames(dat))
  dat2$gene = factor(dat2$gene,levels = rev(colnames(dat)))

  if(!cluster) dat2 = arrange(dat2,group)
  ng = length(unique(group))
  col = legend_color[1:ng]
  names(col) = levels(group)
  p1 = ggplot(dat, aes(x = 1:nrow(dat), y = 1)) +
    geom_tile(aes(fill = group))+
    scale_fill_manual(values = col)+
    theme(panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank()) +
    scale_x_continuous(expand = c(0,0))+
    labs(fill = groupname)

  midpoint = ifelse(fill_mid,stats::median(dat2$exp),0)

  p2 = ggplot(dat2, aes(samples, gene, fill = exp)) +
    geom_tile()+
    theme(panel.grid = element_blank(),
          panel.background = element_rect(fill = NA),
          legend.background = element_rect(fill = NA),
          plot.background = element_rect(fill = NA),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank()
    ) +
    scale_fill_gradient2(low = color[1],
                         mid = color[2],
                         high = color[3],
                         midpoint = midpoint)+
    scale_x_discrete(expand = c(0,0))+
    labs(fill = expname)
  if(!show_rownames) p2 = p2 + theme(axis.text.x = element_blank())
  if(!show_colnames) p2 = p2 + theme(axis.text.y = element_blank())
  p = wrap_plots(p1,p2,nrow = 2,heights = c(1,11))
  return(p)
}

utils::globalVariables(c("gene","samples"))

##' draw_tsne
##'
##' draw tsne plot with annotation by ggplot2
##'
##' @inheritParams draw_pca
##' @param perplexity numeric; perplexity parameter for Rtsne
##' @param color.label color legend label
##' @return a ggplot object
##' @author Xiaojie Sun
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 geom_point
##' @importFrom ggplot2 stat_ellipse
##' @importFrom ggplot2 theme_classic
##' @importFrom ggplot2 theme
##' @importFrom ggplot2 aes
##' @export
##' @examples
##'exp <-  matrix(rnorm(10000),nrow = 50)
##'colnames(exp) <- paste0("sample",1:200)
##'rownames(exp) <- paste0("gene",1:50)
##'exp[1:4,1:4]
##'exp[,1:100] = exp[,1:100]+10
##'group_list <- factor(rep(c("A","B"),each = 100))
##'draw_tsne(exp,group_list)

draw_tsne = function(exp,group_list,perplexity=30,
                     color = c("#2874C5","#f87669","#e6b707","#868686","#92C5DE", "#F4A582", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3",
                               "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3"),
                     color.label = "group",
                     addEllipses = TRUE){
  if(!requireNamespace("Rtsne",quietly = TRUE))
  stop("Package \"Rtsne\" needed for this function to work.
         Please install it by install.packages('Rtsne')",call. = FALSE)
  p1 <-  all(apply(exp,2,is.numeric))
  if(!p1) stop("exp must be a numeric matrix")
  p2  <-  (sum(!duplicated(group_list)) > 1)
  if(!p2) stop("group_list must more than 1")
  p3 <- is.factor(group_list)
  if(!p3) {
    group_list = factor(group_list)
    warning("group_list was covert to factor")
  }
  tsne_out = Rtsne::Rtsne(t(exp),perplexity = perplexity)
  pdat = data.frame(tsne_out$Y,group_list)
  colnames(pdat) = c("Y1","Y2","group")
  p = ggplot(pdat,aes(Y1,Y2))+
    geom_point(aes(Y1,Y2,fill = group),
               shape = 21,color = "black")+
    scale_color_manual(values = color[1:nlevels(group_list)])+
    scale_fill_manual(values = color[1:nlevels(group_list)])+
    theme_classic()+
    theme(legend.position = "top")+
    labs(color = color.label,fill = color.label)
  if(addEllipses) p = p +
    stat_ellipse(aes(color = group,fill = group),
                 geom = "polygon",
                 alpha = 0.3,
                 linetype = 2)
  return(p)
}

utils::globalVariables(c("Y1","Y2","group","patient","ri","time","event"))

##' draw_KM
##'
##' draw KM-plot with two or more group
##'
##' @inheritParams draw_pca
##' @param meta survival data with time and event column
##' @param time_col colname of time
##' @param event_col colname of event
##' @param legend.title legend title
##' @param legend.labs character vector specifying legend labels
##' @param ... other parameters from ggsurvplot
##' @return a KM-plot
##' @author Xiaojie Sun
##' @importFrom survival survfit
##' @importFrom survival Surv
##' @importFrom survminer ggsurvplot
##' @export
##' @examples
##'require("survival")
##'x = survival::lung
##'draw_KM(meta = x,
##'        group_list = x$sex,event_col = "status")

draw_KM = function(meta,
                   group_list,
                   time_col = "time",event_col = "event",
                   legend.title = "Group",
                   legend.labs = levels(group_list),
                   color = c("#2874C5","#f87669","#e6b707","#868686","#92C5DE", "#F4A582", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3",
                             "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3"),
                             ...){
  p1 <-  (time_col %in% colnames(meta)) & (event_col %in% colnames(meta))
  if(!p1){
    stop("meta data must involved time and event columns")
  }else{
    z1 = which(colnames(meta)==time_col)
    z2 = which(colnames(meta)==event_col)
    colnames(meta)[z1] = "time"
    colnames(meta)[z2] = "event"
  }
  p2  <-  (sum(!duplicated(group_list)) > 1)
  if(!p2) stop("group_list must more than 1")
  p3 <- is.factor(group_list)
  if(!p3) {
    group_list = factor(group_list)
    warning("group_list was covert to factor")
  }
  meta$Group = group_list
  sfit <- survival::survfit(Surv(time, event) ~ Group,
                  data = meta)
  p = survminer::ggsurvplot(sfit,pval = T,data = meta,
                 palette = color[1:nlevels(group_list)],
                 legend.title = legend.title,
                 legend.labs = legend.labs,
                 ...)
  return(p$plot)
}

