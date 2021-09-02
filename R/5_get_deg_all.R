##' get_deg_all
##'
##' do diffiencial analysis according to exprission set and group information
##'
##' @inheritParams draw_pca
##' @inheritParams draw_volcano
##' @inheritParams draw_heatmap
##' @inheritParams multi_deg_all
##' @param ids a data.frame with 2 columns,including probe_id and symbol
##' @param entriz logical , if TRUE ,convert symbol to entriz id.
##' @param heat_id id of heatmap,1 for all DEGs,2 for head and tail,3 for top n DEGs
##' @param gene_number how many DEGs will heatmap show .
##' @param color_volcano color for volcano plot
##' @return a list with deg data.frame, volcano plot ,pca plot ,heatmap and a list with DEGs.
##' @author Xiaojie Sun
##' @importFrom patchwork wrap_plots
##' @importFrom ggplot2 ggsave
##' @importFrom stringr str_split
##' @export
##' @examples
##' gse = "GSE42872"
##' geo = geo_download(gse,destdir=tempdir())
##' group_list = rep(c("A","B"),each = 3)
##' group_list = factor(group_list)
##' library(hugene10sttranscriptcluster.db)
##' ids <- toTable(hugene10sttranscriptclusterSYMBOL)
##' dcp = get_deg_all(geo$exp,group_list,ids)
##' head(dcp$deg)
##' dcp$plots
##' @seealso
##' \code{\link{get_deg}};\code{\link{multi_deg_all}}

get_deg_all <- function(exp,
                    group_list,
                    ids,
                    logFC_cutoff=1,
                    pvalue_cutoff=0.05,
                    adjust = FALSE,
                    entriz = TRUE,
                    scale_before = FALSE,
                    n_cutoff = 3,
                    cluster_cols = TRUE,
                    annotation_legend = FALSE,
                    show_rownames = FALSE,
                    legend = FALSE,
                    lab = NA,
                    pkg = 4,
                    symmetry = FALSE,
                    heat_union = TRUE,
                    heat_id = 1,
                    gene_number = 200,
                    color_volcano = c("#2874C5", "grey", "#f87669")) {
  if(nlevels(group_list)==2){
    deg <-  get_deg(exp,group_list,ids,
                    logFC_cutoff=logFC_cutoff,
                    pvalue_cutoff=pvalue_cutoff,
                    adjust = adjust,
                    entriz = entriz)
    cgs = get_cgs(deg)$deg
    volcano_plot = draw_volcano(deg,pkg=pkg,
                                lab =lab,
                                pvalue_cutoff = pvalue_cutoff,
                                logFC_cutoff=logFC_cutoff,
                                adjust = adjust,
                                symmetry = symmetry,
                                color = color_volcano)
    pca_plot = draw_pca(exp,group_list)
    heatmap = draw_heatmap2(exp,group_list,deg,
                            show_rownames = show_rownames,
                            heat_id = heat_id,
                            gene_number = gene_number,
                            scale_before = scale_before,
                            n_cutoff = n_cutoff,
                            cluster_cols = cluster_cols,
                            annotation_legend=annotation_legend
    )
    if(as.numeric(grDevices::dev.cur())!=1) grDevices::graphics.off()
    result = list(
      deg = deg,
      cgs = cgs,
      plots = wrap_plots(heatmap,pca_plot,volcano_plot)+plot_layout(guides = 'collect')
    )
    print(paste0(nrow(cgs$down)," down genes,",nrow(cgs$up)," up genes"))
  }else{
    result <- multi_deg_all(exp,
                            group_list,
                            ids,
                            logFC_cutoff=logFC_cutoff,
                            pvalue_cutoff=pvalue_cutoff,
                            adjust = adjust,
                            entriz = entriz,
                            scale_before = scale_before,
                            n_cutoff = n_cutoff,
                            cluster_cols = cluster_cols,
                            annotation_legend = annotation_legend,
                            show_rownames = show_rownames,
                            legend = legend,
                            lab = lab,
                            pkg = pkg,
                            symmetry = symmetry,
                            heat_union = heat_union,
                            heat_id = heat_id,
                            gene_number = gene_number,
                            color_volcano = color_volcano)
  }
  return(result)
}



