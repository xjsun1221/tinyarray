##' make_tcga_group
##'
##' make tcga group for given tcga expression matrix
##'
##' @inheritParams trans_exp
##' @importFrom stringr str_starts
##' @importFrom stringr str_sub
##' @export
##' @return a group factor with normal and tumor ,conresponse to colnames for expression matrtix
##' @author Xiaojie Sun
##' @examples
##' k = make_tcga_group(exp_hub1);table(k)
##' @seealso
##' \code{\link{geo_download}};\code{\link{draw_volcano}};\code{\link{draw_venn}}

make_tcga_group <- function(exp){
  k1 = stringr::str_starts(colnames(exp),"TCGA")
  if(!any(k1))stop("no tcga samples detected,please check it")
  k2 = as.numeric(stringr::str_sub(colnames(exp),14,15))<10
  group_list = ifelse(k1&k2,"tumor","normal")
  group_list = factor(group_list,levels = c("normal","tumor"))
  return(group_list)
}

##' exp_surv
##'
##' draw surv plot for a hub gene expression matrix for tumor samples
##'
##' @inheritParams surv_KM
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
##' \code{\link{geo_download}};\code{\link{draw_volcano}};\code{\link{draw_venn}}

exp_surv <- function(exprSet_hub,meta){
  cut.point = point_cut(exprSet_hub,meta)
  splots <- lapply(rownames(exprSet_hub), function(g){
    i = which(rownames(exprSet_hub)== g)
    meta$gene=ifelse(as.integer(exprSet_hub[g,]) > cut.point[[i]],'high','low')
    sfit1=survival::survfit(survival::Surv(time, event)~gene, data=meta)
    p = survminer::ggsurvplot(sfit1,pval =TRUE,
                              palette = c("red", "grey"),
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
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 aes_string
##' @importFrom ggplot2 geom_boxplot
##' @importFrom ggplot2 scale_fill_manual
##' @importFrom ggpubr stat_compare_means
##' @importFrom ggplot2 theme_classic
##' @importFrom ggplot2 theme
##' @export
##' @return box plots list for all genes in the matrix
##' @author Xiaojie Sun
##' @examples
##' k = exp_boxplot(log2(exp_hub1+1));k[[1]]
##' @seealso
##' \code{\link{geo_download}};\code{\link{draw_volcano}};\code{\link{draw_venn}}

exp_boxplot <-  function(exp_hub){
  if(nrow(exp_hub)>15)warning(paste0(nrow(exp_hub)," figures will be produced"))
  group_list = make_tcga_group(exp_hub)
  dat = as.data.frame(t(exp_hub))
  dat$group_list = group_list
  bplots = lapply(1:nrow(exp_hub),function(i){
    ggplot(dat,
           aes_string(x = "group_list",
                      y = rownames(exp_hub)[[i]],
                      fill = "group_list"))+
      geom_boxplot()+
      scale_fill_manual(values = c("normal" = "grey",
                                   "tumor" = "red"))+
      stat_compare_means(method = "t.test")+
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
##' \code{\link{geo_download}};\code{\link{draw_volcano}};\code{\link{draw_venn}}

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
