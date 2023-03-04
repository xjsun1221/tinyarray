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
##' @param color_volcano color for volcano plot
##' @param ... other parameters from get_deg
##' @return a list with deg data.frame, volcano plot ,pca plot ,heatmap and a list with DEGs.
##' @author Xiaojie Sun
##' @importFrom patchwork wrap_plots
##' @importFrom ggplot2 ggsave
##' @importFrom stringr str_split
##' @export
##' @examples
##' \donttest{
##' gse = "GSE42872"
##' geo = geo_download(gse,destdir=tempdir())
##' group_list = rep(c("A","B"),each = 3)
##' group_list = factor(group_list)
##' find_anno(geo$gpl)
##' ids <- AnnoProbe::idmap(geo$gpl,destdir = tempdir())
##' dcp = get_deg_all(geo$exp,group_list,ids)
##' head(dcp$deg)
##' dcp$plots
##' }
##' @seealso
##' \code{\link{get_deg}};\code{\link{multi_deg_all}}

get_deg_all <- function(exp,
                    group_list,
                    ids,pkg=4,
                    color_volcano = c("#2874C5", "grey", "#f87669"),
                    my_genes = NULL,
                    show_rownames = FALSE,
                    entriz = TRUE,
                    ...) {
  if(nlevels(group_list)==2){
    deg <-  get_deg(exp,group_list,ids,entriz = entriz,...)
    cgs = get_cgs(deg)$deg
    volcano_plot = draw_volcano2(deg,pkg = pkg,color = color_volcano)
    pca_plot = draw_pca(exp,group_list)
    heatmap = draw_heatmap2(exp,group_list,deg,my_genes,show_rownames = show_rownames)
    if(as.numeric(grDevices::dev.cur())!=1) grDevices::graphics.off()
    result = list(
      deg = deg,
      cgs = cgs,
      plots = wrap_plots(heatmap,pca_plot,volcano_plot)+plot_layout(guides = 'collect')
    )
    message(paste0(nrow(cgs$down)," down genes,",nrow(cgs$up)," up genes"))
  }else{
    result <- multi_deg_all(exp,
                            group_list,
                            ids,
                            pkg = pkg,
                            ...)
  }
  return(result)
}



