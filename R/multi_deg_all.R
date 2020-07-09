##' multi_deg_all
##'
##' do diffiencial analysis according to exprission set and group information
##'
##' @inheritParams draw_pca
##' @inheritParams draw_volcano
##' @inheritParams draw_heatmap
##' @inheritParams draw_heatmap2
##' @inheritParams multi_deg
##' @return a list with deg data.frame, volcano plot and a list with DEGs.
##' @author Xiaojie Sun
##' @importFrom patchwork wrap_plots
##' @importFrom stringr str_split
##' @importFrom dplyr union_all
##' @importFrom patchwork plot_layout
##' @export
##' @examples
##' gse = "GSE474"
##' geo = geo_download(gse)
##' geo$exp[1:4,1:4]
##' geo$exp=log2(geo$exp+1)
##' library(stringr)
##' group_list=ifelse(str_detect(geo$pd$title,"MObese"),"MObese",ifelse(str_detect(geo$pd$title,"NonObese"),"NonObese","Obese"))
##' group_list=factor(group_list,levels = c("NonObese","Obese","MObese"))
##' find_anno(geo$gpl)
##' ids <- toTable(hgu133aSYMBOL)
##' dcp = multi_deg_all(geo$exp,group_list,ids,adjust = F,heat_id = 2,gene_number = 360)
##' dcp[[3]]
##' @seealso
##' \code{\link{geo_download}};\code{\link{draw_volcano}};\code{\link{draw_venn}}

multi_deg_all <- function(exp,
                          group_list,
                          ids,
                          logFC_cutoff=1,
                          pvalue_cutoff=0.05,
                          adjust = F,
                          entriz = T,
                          scale_before = F,
                          n_cutoff = 2,
                          cluster_cols = T,
                          annotation_legend = F,
                          lab = NA,
                          pkg = 4,
                          symmetry = F,
                          heat_union = T,
                          heat_id = 2,
                          gene_number = 200) {
  deg = multi_deg(
    exp,
    group_list,
    ids,
    logFC_cutoff = logFC_cutoff,
    pvalue_cutoff = pvalue_cutoff,
    adjust = adjust,
    entriz = entriz
  )
  #exp = data.frame(exp)
  #exp = exp[match(deg[[1]]$probe_id,rownames(exp)),]
  cgs = get_cgs(deg)
  volcano_plot = draw_volcano2(deg)
  pca_plot = draw_pca(exp,group_list)
  heatmap = draw_heatmap2(geo$exp,group_list,
                          deg,
                          heat_id=heat_id,
                          gene_number=gene_number,
                          scale_before = scale_before,
                          n_cutoff = n_cutoff,
                          cluster_cols = cluster_cols,
                          annotation_legend=annotation_legend
                          )
  x = lapply(cgs,function(x)x$diff$diffprobes)
  venn = draw_venn(x," ")
  if(as.numeric(dev.cur())!=1) graphics.off()
  plotlist = list(heatmap,pca_plot,venn,volcano_plot)
  layout <- '
  AABBCC
  AABBCC
  DDDDDD
  DDDDDD
  '
  result = list(
    deg = deg,
    cgs = cgs,
    plots = wrap_plots(plotlist) +
      plot_layout(design = layout) +
      plot_layout(guides = 'collect')
  )
  diffprobes = lapply(cgs,function(x)x$diff$diffprobes)
  print(paste0(length(union_all(diffprobes))," DEGs in all,",length(intersect_all(diffprobes))," DEGs in common."))
  return(result)
  }
