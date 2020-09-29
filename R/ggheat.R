##' ggheat
##'
##' draw heatmap plot with annotation by ggplot2
##'
##' @param dat expression matrix for plot
##' @param group group for expression colnames
##' @param cluster logical,cluster or not, default F
##' @param show_rownames logical,show rownames in plot or not,default T
##' @param show_colnames logical,show colnames in plot or not,default T
##' @param groupname name of group legend
##' @param expname name of exp legene
##' @return a ggplot object
##' @author Xiaojie Sun
##' @importFrom pheatmap pheatmap
##' @importFrom tibble rownames_to_column
##' @importFrom tidyr gather
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
##' rm(list = ls())
##' exp_dat = matrix(sample(100:1000,40),ncol = 4)
##' exp_dat[1:(nrow(exp_dat)/2),] =  exp_dat[1:(nrow(exp_dat)/2),]-1000
##' rownames(exp_dat) = paste0("sample",1:nrow(exp_dat))
##' colnames(exp_dat) = paste0("gene",1:ncol(exp_dat))

##' group = rep(c("A","B"),each = nrow(exp_dat)/2)
##' group = factor(group,levels = c("A","B"))
##' ggheat(exp_dat,group)
##' ggheat(exp_dat,group,cluster = T)
##' ggheat(exp_dat,group,cluster = T,show_rownames = F,
##'        show_colnames = F,groupname = "risk",expname = "expression")



ggheat = function(dat,group,cluster = F,
                  show_rownames = T,show_colnames = T,
                  groupname = "group",expname = "exp"){
  dat = data.frame(dat)

  if(cluster){
    ph = pheatmap::pheatmap(t(dat),silent = T)
    dat = dat[ph$tree_col$order,ph$tree_row$order]
    group = group[ph$tree_col$order]
  }

  dat$group = group
  dat2 = gather(rownames_to_column(dat,
                                       "samples"),
                    gene,exp,-group,-samples)

  dat2$samples = factor(dat2$samples,levels = rownames(dat))
  dat2$gene = factor(dat2$gene,levels = rev(colnames(dat)))

  if(!cluster) dat2 = arrange(dat2,group)
  col = c("blue","red")
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
    scale_fill_gradient2(low = "blue",
                         mid = "white",
                         high = "red",
                         midpoint = median(dat2$exp))+
    scale_x_discrete(expand = c(0,0))+
    labs(fill = expname)
  if(!show_rownames) p2 = p2 + theme(axis.text.x = element_blank())
  if(!show_colnames) p2 = p2 + theme(axis.text.y = element_blank())
  p = wrap_plots(p1,p2,nrow = 2,heights = c(1,11))
  return(p)
}
