##' exp_surv
##'
##' draw surv plot for a hub gene expression matrix for tumor samples
##'
##' @inheritParams surv_KM
##' @inheritParams exp_boxplot
##' @importFrom survival survfit
##' @importFrom survival Surv
##' @importFrom survminer ggsurvplot
##' @importFrom ggplot2 theme
##' @importFrom ggplot2 element_text
##' @export
##' @return survival plots list for all genes
##' @author Xiaojie Sun
##' @examples
##' tmp = exp_surv(exprSet_hub1,meta1)
##' patchwork::wrap_plots(tmp)
##' @seealso
##' \code{\link{exp_boxplot}};\code{\link{box_surv}};\code{\link{draw_venn}}

exp_surv <- function(exprSet_hub,meta,color = c("grey", "red")){
  cut.point = point_cut(exprSet_hub,meta)
  splots <- lapply(rownames(exprSet_hub), function(g){
    i = which(rownames(exprSet_hub)== g)
    meta$gene=ifelse(as.numeric(exprSet_hub[g,]) > cut.point[[i]],'high','low')
    if(length(unique(meta$gene))==1) stop(paste0("gene",g,"with too low expression"))
    sfit1=survival::survfit(survival::Surv(time, event)~gene, data=meta)
    p = survminer::ggsurvplot(sfit1,pval =TRUE,
                              palette = rev(color),
                              data = meta,
                              legend = c(0.8,0.8),
                              title = rownames(exprSet_hub)[[i]]
    )
    p2 = p$plot+
      theme(plot.title = element_text(hjust = 0.5))
    return(p2)
  })
  return(splots)
}

##' exp_boxplot
##'
##' draw box plot for a hub gene expression matrix
##'
##' @param exp_hub an expression matrix for hubgenes
##' @param color color for boxplot
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 aes_string
##' @importFrom ggplot2 geom_boxplot
##' @importFrom ggplot2 scale_fill_manual
##' @importFrom ggplot2 theme_classic
##' @importFrom ggplot2 theme
##' @export
##' @return box plots list for all genes in the matrix
##' @author Xiaojie Sun
##' @examples
##' k = exp_boxplot(log2(exp_hub1+1));k[[1]]
##' @seealso
##' \code{\link{exp_surv}};\code{\link{box_surv}}

exp_boxplot <-  function(exp_hub,color = c("grey","red")){
  if(nrow(exp_hub)>15)warning(paste0(nrow(exp_hub)," figures will be produced"))
  group_list = make_tcga_group(exp_hub)
  k = rownames(exp_hub)
  rownames(exp_hub) = paste0("gene",1:nrow(exp_hub))
  dat = as.data.frame(t(exp_hub))
  dat$group_list = group_list
  color = color[1:length(levels(group_list))]
  names(color) = levels(group_list)
  bplots = lapply(1:nrow(exp_hub),function(i){
    ggplot(dat,
           aes_string(x = "group_list",
                      y = rownames(exp_hub)[[i]],
                      fill = "group_list"))+
      geom_boxplot()+
      scale_fill_manual(values = color)+
      labs(y = k[[i]])+
      ggpubr::stat_compare_means(method = "t.test")+
      #stat_compare_means(label.y = 15)  +
      theme_classic()+
      theme(legend.position = "none")
  })
  return(bplots)
}

##' box_surv
##'
##' draw box plot for a hub gene expression matrix
##'
##' @inheritParams exp_boxplot
##' @inheritParams exp_surv
##' @importFrom patchwork plot_layout
##' @importFrom ggplot2 theme
##' @export
##' @return patchwork result for hub genes boxplot and survival plot
##' @author Xiaojie Sun
##' @examples
##' k = box_surv(log2(exp_hub1+1),exprSet_hub1,meta1);k[[1]]
##' @seealso
##' \code{\link{exp_boxplot}};\code{\link{exp_surv}}

box_surv <-function(exp_hub,exprSet_hub,meta){
  splots = exp_surv(exprSet_hub,meta)
  boxplots = exp_boxplot(exp_hub)
  layout <- '
  ABB
  '
  wp = lapply(1:nrow(exprSet_hub),function(i){
    boxplots[[i]]+splots[[i]]+
      plot_layout(design = layout)& theme(legend.position = "none")
  })
  return(wp)
}
