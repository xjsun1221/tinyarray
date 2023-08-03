##' plot_deg
##'
##' plot pca plot,volcano plot,heatmap,and venn plot for  Differential analysis result
##'
##' @inheritParams multi_deg_all
##' @return plots
##' @author Xiaojie Sun
##' @export
##' @examples
##' \dontrun{
##' gse = "GSE474"
##' geo = geo_download(gse,destdir=tempdir())
##' geo$exp[1:4,1:4]
##' geo$exp=log2(geo$exp+1)
##' group_list=ifelse(stringr::str_detect(geo$pd$title,"MObese"),"MObese",
##' ifelse(stringr::str_detect(geo$pd$title,"NonObese"),"NonObese","Obese"))
##' group_list=factor(group_list,levels = c("NonObese","Obese","MObese"))
##' find_anno(geo$gpl)
##' ids = AnnoProbe::idmap(geo$gpl,destdir = tempdir())
##' deg = get_deg(geo$exp,group_list,ids,adjust = FALSE)
##' plot_deg(geo$exp,group_list,deg)
##' }
plot_deg = function(exp,
                     group_list,
                     symmetry = TRUE,
                     my_genes = NULL,
                     show_rownames = FALSE,
                     cluster_cols = TRUE,
                     color_volcano = c("#2874C5", "grey", "#f87669"),
                     pvalue_cutoff = 0.05,
                     logFC_cutoff = 1,
                     annotation_legend = FALSE,
                     lab = NA,
                     species = "human"
){
  cgs = get_cgs(deg)
  volcano_plot = draw_volcano2(deg,
                               pkg = 4,
                               symmetry = symmetry,
                               color = color_volcano,
                               pvalue_cutoff = pvalue_cutoff,
                               logFC_cutoff = logFC_cutoff,
                               adjust = adjust
  )
  pca_plot = draw_pca(exp,group_list)
  heatmap = draw_heatmap2(exp,group_list,deg,my_genes,
                          show_rownames = show_rownames,
                          cluster_cols = cluster_cols)
  x = lapply(cgs,function(x)x$diff$diffprobes)
  venn = draw_venn(x," ")
  if(as.numeric(grDevices::dev.cur())!=1) grDevices::graphics.off()
  plotlist = list(heatmap,pca_plot,venn,volcano_plot)
  layout <- '
  AABBCC
  AABBCC
  DDDDDD
  DDDDDD
  '
  plots = wrap_plots(plotlist) +
    plot_layout(design = layout) +
    plot_layout(guides = 'collect')
  return(plots)
}





