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
##' patchwork::wrap_plots(tmp)+patchwork::plot_layout(guides = "collect")
##' tmp2 = exp_surv(exprSet_hub1,meta1,cut.point = TRUE)
##' patchwork::wrap_plots(tmp2)+patchwork::plot_layout(guides = "collect")
##' @seealso
##' \code{\link{exp_boxplot}};\code{\link{box_surv}};\code{\link{draw_venn}}


exp_surv <- function(exprSet_hub,meta,cut.point = FALSE,color = c("#2874C5", "#f87669")){
  if(cut.point)cut_point = point_cut(exprSet_hub,meta)

  splots <- lapply(rownames(exprSet_hub), function(g){
    i = which(rownames(exprSet_hub)== g)
    if(cut.point){
      meta$gene=ifelse(as.numeric(exprSet_hub[g,]) > cut_point[[i]],'high','low')
    }else{
      meta$gene=ifelse(as.numeric(exprSet_hub[g,]) > stats::median(as.numeric(exprSet_hub[g,])),'high','low')
    }
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



##' risk_plot
##'
##' draw risk plot
##'
##' @inheritParams exp_surv
##' @param riskscore a numeric vector of riskscore
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 aes
##' @importFrom ggplot2 geom_point
##' @importFrom ggplot2 scale_color_manual
##' @importFrom ggplot2 geom_vline
##' @importFrom ggplot2 scale_x_continuous
##' @importFrom ggplot2 theme_bw
##' @importFrom ggplot2 theme
##' @importFrom ggplot2 element_blank
##' @importFrom ggplot2 labs
##' @importFrom ggplot2 scale_color_manual
##' @importFrom ggplot2 scale_y_continuous
##' @importFrom ggplot2 element_line
##' @importFrom patchwork wrap_plots
##' @export
##' @return risk plot
##' @author Xiaojie Sun
##' @examples
##' risk_plot(exprSet_hub1,meta1,riskscore = rnorm(nrow(meta1)))
##' @seealso
##' \code{\link{exp_boxplot}};\code{\link{box_surv}};\code{\link{draw_venn}}

risk_plot = function(exprSet_hub,meta,riskscore,
                     cut.point = FALSE,color = c("#2fa1dd","#f87669")){
  if(ncol(exprSet_hub) != nrow(meta))stop("your exprSet_hub is not corresponds to meta")
  if(length(riskscore) != nrow(meta))stop("riskscore is not corresponds to meta")
  if (nrow(exprSet_hub)>30) {
    warning("seems too many of genes in heatmap")
  }
  if (nrow(exprSet_hub)>100) {
    stop("too many of genes in heatmap")
  }
  meta$riskscorefp = riskscore
  if(cut.point){
    cut = survminer::surv_cutpoint(
      meta,
      time = "time",
      event = "event",
      variables = "riskscorefp")[["cutpoint"]][1,1]
  }else{
    cut = stats::median(riskscore)
  }
  meta$Risk = ifelse(riskscore>cut,'high','low')
  meta$Risk = factor(meta$Risk,levels = c("low","high"))
  fp_dat=data.frame(patientid=1:length(riskscore),
                    riskscore=as.numeric(sort(riskscore)),
                    ri = meta$Risk[order(riskscore)])
  sur_dat=data.frame(patientid=1:length(riskscore),
                     time=meta[order(riskscore),'time'] ,
                     event=meta[order(riskscore),'event'])
  sur_dat$event=ifelse(sur_dat$event==0,'alive','death')

  # exp

  exp_dat=scale(t(exprSet_hub[,order(riskscore)]))

  ###第一个图----
  p1=ggplot(fp_dat,aes(x=patientid,y=riskscore,color = ri))+
    geom_point()+
    scale_color_manual(values = color)+
    geom_vline(xintercept = sum(riskscore<cut),lty = 2)+
    scale_x_continuous(expand=c(0,0))+
    theme_bw()+
    theme(legend.position = "none",
          axis.line = element_line(color='black'),
          plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    labs(x = "",y = "Riskscore")
  p1
  #第二个图
  p2=ggplot(sur_dat,aes(x=patientid,y=time))+
    geom_point(aes(col=event))+
    geom_vline(xintercept = sum(riskscore<cut),lty = 2)+
    scale_color_manual(values = color)+
    scale_x_continuous(expand=c(0,0))+
    scale_y_continuous(expand=c(0,0))+
    theme_bw()+
    labs(color = "Event",y = "Time",x = "")+
    theme(axis.line = element_line(color='black'),
          plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  p2
  #第三个图
  p3 = ggheat(exp_dat,
              meta$Risk[order(riskscore)],
              show_rownames = F,
              color = c(color[[1]],"white",color[[2]]),
              legend_color = color,groupname = "Risk",expname = "Expression")

  n = ifelse(nrow(exprSet_hub)<15,3,
             ifelse(nrow(exprSet_hub)<20,4,5))
  wrap_plots(p1,p2,p3,heights = c(1,1,n))
}
utils::globalVariables(c("patientid","genekitr"))
