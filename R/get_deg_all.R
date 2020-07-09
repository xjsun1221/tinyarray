##' get_deg_all
##'
##' do diffiencial analysis according to exprission set and group information
##'
##' @inheritParams draw_pca
##' @inheritParams draw_volcano
##' @inheritParams draw_heatmap
##' @param ids a data.frame with 2 columns,including probe_id and symbol
##' @param entriz logical , if TRUE ,convert symbol to entriz id.
##' @param heat_id id of heatmap,1 for all DEGs,2 for head and tail,3 for top n DEGs
##' @param gene_number how many DEGs will heatmap show .
##' @return a list with deg data.frame, volcano plot and a list with DEGs.
##' @author Xiaojie Sun
##' @importFrom patchwork wrap_plots
##' @importFrom ggplot2 ggsave
##' @importFrom stringr str_split
##' @export
##' @examples
##' gse = "GSE42872"
##' geo = geo_download(gse)
##' group_list = rep(c("A","B"),each = 3)
##' ids = AnnoProbe::idmap('GPL6244')
##' get_deg_all(geo$exp,group_list,ids)
##' @seealso
##' \code{\link{geo_download}};\code{\link{draw_volcano}};\code{\link{draw_venn}}

get_deg_all <- function(exp,
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
                    heat_id = 2,
                    gene_number = 200) {
  deg <-  get_deg(exp,group_list,ids,
                logFC_cutoff=logFC_cutoff,
                pvalue_cutoff=pvalue_cutoff,
                adjust = adjust,
                entriz = entriz)
  cgs <- list(upgenes = deg$symbol[deg$change=="up"],
              downgenes = deg$symbol[deg$change=="down"],
              diffgenes = deg$symbol[deg$change!="stable"])
  cd = gene_number> length(cgs$upgenes) | gene_number> length(cgs$downgenes)
  er = paste0("gene_number must less than ",min(c(length(cgs$upgenes),length(cgs$downgenes))))
  if(heat_id==3 & cd) stop(er)
  if(heat_id==2 & gene_number> length(cgs$diffgenes))stop(
    paste0("gene_number must less than ",length(cgs$diffgenes))
  )
  volcano_plot = draw_volcano(deg,pkg=pkg,
                              lab =lab,
                              pvalue_cutoff = pvalue_cutoff,
                              logFC_cutoff=logFC_cutoff,
                              adjust = adjust,
                              symmetry = symmetry)
  pca_plot = draw_pca(exp,group_list)
  heatmap = draw_heatmap2(exp,group_list,deg,
                          heat_id = heat_id,
                          gene_number = gene_number,
                          scale_before = scale_before,
                          n_cutoff = n_cutoff,
                          cluster_cols = cluster_cols,
                          annotation_legend=annotation_legend
                          )
  if(as.numeric(dev.cur())!=1) graphics.off()
  result = list(
    deg = deg,
    cgs = cgs,
    plots = wrap_plots(heatmap,pca_plot,volcano_plot)+plot_layout(guides = 'collect')
  )
  print(paste0(length(cgs$downgenes)," down genes,",length(cgs$upgenes)," up genes"))
  return(result)
}

##' find annotation package or files
##'
##' find gpl annotation package or files
##'
##' @param gpl a gpl accession
##' @param install if R packages will be installed
##' @return a list with deg data.frame, volcano plot and a list with DEGs.
##' @author Xiaojie Sun
##' @importFrom stringr str_remove_all
##' @importFrom stringr str_to_upper
##' @export
##' @examples
##' find_anno(GPL570)
##' @seealso
##' \code{\link{geo_download}};\code{\link{draw_volcano}};\code{\link{draw_venn}}

find_anno <-function(gpl,install = F){
  gpl = str_to_upper(gpl)
  data("pkg_all")
  data("exists_anno_list")
  if(!any(pkg_all$gpl==gpl)) {
    # R包不可用
    if(gpl %in% setdiff(exists_anno_list,pkg_all$gpl)){
      # 只有idmap可用
      ml1 = str_remove_all(paste0("`ids <- AnnoProbe::idmap\\(","\\'",gpl,"\\'","\\)`"),"\\\\")
      print(paste0("no annotation packages avliable,please use ",ml1))
    }else{
      # R包和idmap都不可用
      print("no annotation avliable in Bioconductor and AnnoProbe")
    }
  }else {
    qz = pkg_all$bioc_package[pkg_all$gpl== gpl]
    ml1 = str_remove_all(paste0("`ids <- AnnoProbe::idmap\\(","\\'",gpl,"\\'","\\)`"),"\\\\")
    ml2 = str_remove_all(paste0("`ids <- toTable\\(",qz,"SYMBOL\\)`"),"\\\\")
    if(install){
      if(!suppressMessages(require(paste0(qz,".db"),character.only = T)))BiocManager::install(paste0(qz,".db"))
      suppressMessages(library(paste0(qz,".db"),character.only = T))
    }
    if(!(gpl %in% exists_anno_list)) {
      #仅有R包可用
      print(paste0(ml2," is avaliable"))
    }else {
      #idmap和R包都可用
      print(paste0(ml2," and ",ml1 ," are both avaliable"))
      }
  }
}


