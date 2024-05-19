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
##' @return a list with deg data.frame, volcano plot ,pca plot ,heatmap and a list with DEGs.
##' @author Xiaojie Sun
##' @importFrom patchwork wrap_plots
##' @importFrom ggplot2 ggsave
##' @importFrom stringr str_split
##' @export
##' @examples
##' \dontrun{
##' if(requireNamespace("Biobase",quietly = TRUE)&
##'    requireNamespace("AnnoProbe",quietly = TRUE)){
##'   gse = "GSE42872"
##'   geo = geo_download(gse,destdir=tempdir())
##'   group_list = rep(c("A","B"),each = 3)
##'   group_list = factor(group_list)
##'   find_anno(geo$gpl)
##'   ids <- AnnoProbe::idmap(geo$gpl,destdir = tempdir())
##'   dcp = get_deg_all(geo$exp,group_list,ids,entriz = FALSE)
##'   head(dcp$deg)
##'   dcp$plots
##' }else{
##'   if(!requireNamespace("AnnoProbe",quietly = TRUE)) {
##'     warning("Package 'AnnoProbe' needed for this function to work.
##'          Please install it by install.packages('AnnoProbe')",call. = FALSE)
##'   }
##'   if(!requireNamespace("Biobase",quietly = TRUE)) {
##'     warning("Package 'Biobase' needed for this function to work.
##'          Please install it by BiocManager::install('Biobase')",call. = FALSE)
##'   }
##' }
##' }
##' @seealso
##' \code{\link{get_deg}};\code{\link{multi_deg_all}}

get_deg_all <- function(exp,
                        group_list,
                        ids,
                        symmetry = TRUE,
                        my_genes = NULL,
                        show_rownames = FALSE,
                        cluster_cols = TRUE,
                        color_volcano = c("#2874C5", "grey", "#f87669"),
                        logFC_cutoff=1,
                        pvalue_cutoff=0.05,
                        adjust = FALSE,
                        entriz = TRUE,
                        n_cutoff = 2,
                        annotation_legend = FALSE,
                        lab = NA,
                        species = "human") {
  if(nlevels(group_list)==2){
    deg <-  get_deg(exp,group_list,ids,
                    logFC_cutoff=logFC_cutoff,
                    pvalue_cutoff=pvalue_cutoff,
                    adjust = adjust,
                    entriz = entriz,
                    species = species)
    cgs = get_cgs(deg)
    volcano_plot = draw_volcano(deg,pkg=4,
                                lab =lab,
                                pvalue_cutoff = pvalue_cutoff,
                                logFC_cutoff=logFC_cutoff,
                                adjust = adjust,
                                symmetry = symmetry)
    pca_plot = draw_pca(exp,group_list)
    heatmap = draw_heatmap2(exp,group_list,deg,my_genes,
                            show_rownames = show_rownames,
                            n_cutoff = n_cutoff,
                            cluster_cols = cluster_cols,
                            annotation_legend=annotation_legend)
    if(as.numeric(grDevices::dev.cur())!=1) grDevices::graphics.off()
    result = list(
      deg = deg,
      cgs = cgs,
      plots = wrap_plots(heatmap,pca_plot,volcano_plot)+plot_layout(guides = 'collect')
    )
    message(paste0(nrow(cgs$deg$down)," down genes,",nrow(cgs$deg$up)," up genes"))
  }else{
    result <- multi_deg_all(exp,
                            group_list,
                            ids,
                            logFC_cutoff = logFC_cutoff,
                            pvalue_cutoff = pvalue_cutoff,
                            symmetry = symmetry,
                            my_genes = my_genes,
                            show_rownames = show_rownames,
                            cluster_cols = cluster_cols,
                            color_volcano = color_volcano,
                            adjust = adjust,
                            entriz = entriz,
                            species = species)
  }
  return(result)
}



