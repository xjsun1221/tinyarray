##' get_cgs
##'
##' extract DEGs from deg data.frame
##'
##' @inheritParams draw_volcano
##' @return a list with upgenes,downgenes,diffgenes.
##' @author Xiaojie Sun
##' @export
##' @examples
##' \donttest{
##' #two group
##' gse = "GSE42872"
##' geo = geo_download(gse,destdir=tempdir(),by_annopbrobe = FALSE)
##' group_list = rep(c("A","B"),each = 3)
##' ids = AnnoProbe::idmap('GPL6244',destdir=tempdir())
##' deg = get_deg(geo$exp,group_list,ids)
##' cgs = get_cgs(deg)
##' #mutigroup
##' gse = "GSE474"
##' geo = geo_download(gse,destdir=tempdir(),by_annopbrobe = FALSE)
##' geo$exp[1:4,1:4]
##' geo$exp=log2(geo$exp+1)
##' group_list=ifelse(stringr::str_detect(geo$pd$title,"MObese"),"MObese",
##' ifelse(stringr::str_detect(geo$pd$title,"NonObese"),"NonObese","Obese"))
##' group_list=factor(group_list,levels = c("NonObese","Obese","MObese"))
##' find_anno(geo$gpl)
##' ids = AnnoProbe::idmap(geo$gpl,destdir = tempdir())
##' deg = multi_deg(geo$exp,group_list,ids,adjust = FALSE)
##' cgs = get_cgs(deg)
##' }
##' @seealso
##' \code{\link{geo_download}};\code{\link{draw_volcano}};\code{\link{draw_venn}}


get_cgs <- function(deg){
  if(!is.list(deg) & is.data.frame(deg))stop("deg is a data.frame or list returned by limma")
  cgs = list()
  if(is.data.frame(deg)) deg = list(deg = deg)
  for(i in 1:length(deg)){
    cgs[[i]] = list(up = data.frame(upgenes =deg[[i]]$symbol[deg[[i]]$change=="up"],
                                    upprobes = deg[[i]]$probe_id[deg[[i]]$change=="up"],
                                    stringsAsFactors = FALSE),
                    down = data.frame(downgenes = deg[[i]]$symbol[deg[[i]]$change=="down"],
                                      downprobes = deg[[i]]$probe_id[deg[[i]]$change=="down"],
                                      stringsAsFactors = FALSE),
                    diff = data.frame(diffgenes = deg[[i]]$symbol[deg[[i]]$change!="stable"],
                                      diffprobes = deg[[i]]$probe_id[deg[[i]]$change!="stable"],
                                      stringsAsFactors = FALSE)
    )
  }
  if(!is.data.frame(deg)) names(cgs) = names(deg)
  if(is.data.frame(deg)) cgs = cgs[[1]]
  return(cgs)
}

##' draw_volcano2
##'
##' print one or more volcano plot for Differential analysis result in data.frame fomat.
##'
##' @inheritParams draw_volcano
##' @return one or more volcano plot
##' @author Xiaojie Sun
##' @importFrom patchwork wrap_plots
##' @importFrom patchwork plot_layout
##' @export
##' @examples
##' \donttest{
##' #two group
##' gse = "GSE42872"
##' geo = geo_download(gse,destdir=tempdir(),by_annopbrobe = FALSE)
##' group_list = rep(c("A","B"),each = 3)
##' ids = AnnoProbe::idmap('GPL6244',destdir = tempdir())
##' deg = get_deg(geo$exp,group_list,ids)
##' draw_volcano2(deg)
##' #multigroup
##' gse = "GSE474"
##' geo = geo_download(gse,destdir=tempdir(),by_annopbrobe = FALSE)
##' geo$exp[1:4,1:4]
##' geo$exp=log2(geo$exp+1)
##' group_list=ifelse(stringr::str_detect(geo$pd$title,"MObese"),"MObese",
##' ifelse(stringr::str_detect(geo$pd$title,"NonObese"),"NonObese","Obese"))
##' group_list=factor(group_list,levels = c("NonObese","Obese","MObese"))
##' find_anno(geo$gpl)
##' ids <- AnnoProbe::idmap(geo$gpl,destdir = tempdir())
##' deg = multi_deg(geo$exp,group_list,ids,adjust = FALSE)
##' draw_volcano2(deg)
##' draw_volcano2(deg,color = c("darkgreen","grey","darkred"))
##' }
##' @seealso
##' \code{\link{geo_download}};\code{\link{draw_volcano}};\code{\link{draw_venn}}

draw_volcano2 = function(deg,
                         pkg=4,
                         lab,
                         ...
){
  if(!is.list(deg) & is.data.frame(deg))stop("deg is a data.frame or list returned by limma")
  if(is.data.frame(deg)) {
    volcano_plots = draw_volcano(
      deg ,
      pkg = pkg,
      lab = "logFC",
      ...
    )
  }else{
    volcano_plots <- lapply(1:length(deg),
                            function(k) {
                              draw_volcano(
                                deg[[k]] ,
                                pkg = pkg,
                                lab = names(deg)[k],
                                ...
                              )
                            })
    volcano_plots = wrap_plots(volcano_plots)+
      plot_layout(design = paste(rep(LETTERS[1:length(deg)]),collapse = "")) +
      plot_layout(guides = 'collect')
  }
  return(volcano_plots)
}

##' draw heatmap plots
##'
##' print heatmap plots for expression matrix and group by group_list paramter
##'
##' @inheritParams draw_volcano
##' @inheritParams draw_heatmap
##' @inheritParams draw_pca
##' @param heat_union logical ,use union or intersect DEGs for heatmap
##' @param my_genes genes for pheatmap
##' @return a heatmap plot according to \code{exp} and grouped by \code{group}.
##' @author Xiaojie Sun
##' @importFrom pheatmap pheatmap
##' @export
##' @examples
##' \donttest{
##' gse = "GSE474"
##' geo = geo_download(gse,destdir=tempdir(),by_annopbrobe = FALSE)
##' geo$exp[1:4,1:4]
##' geo$exp=log2(geo$exp+1)
##' group_list=ifelse(stringr::str_detect(geo$pd$title,"MObese"),"MObese",
##' ifelse(stringr::str_detect(geo$pd$title,"NonObese"),"NonObese","Obese"))
##' group_list=factor(group_list,levels = c("NonObese","Obese","MObese"))
##' find_anno(geo$gpl)
##' ids <- AnnoProbe::idmap(geo$gpl,destdir = tempdir())
##' deg = multi_deg(geo$exp,group_list,ids,adjust = FALSE)
##' draw_heatmap2(geo$exp,group_list,deg)
##' }
##' @seealso
##' \code{\link{draw_pca}};\code{\link{draw_volcano}};\code{\link{draw_venn}}

draw_heatmap2 <- function(exp,
                          group_list,
                          deg,
                          my_genes = NULL,
                          heat_union = TRUE,...
){
  if(!is.data.frame(deg)){
    m = deg[[1]]
  }else{
    m = deg
  }
  if(all(m$probe_id %in% rownames(exp))){
    exp = exp[m$probe_id,]
    rownames(exp) = m$symbol
  }else if(all(m$symbol %in% rownames(exp))) {
    exp = exp[m$symbol,]
  }
  if(is.null(my_genes)){
    cgs = get_cgs(deg)
    if(length(cgs)==1){
      cg = cgs[[1]]$diff$diffgenes
    }else{
      if(heat_union){
        cg = union_all(lapply(cgs,function(x)x$diff$diffgenes))
      }else{
        cg = intersect_all(lapply(cgs,function(x)x$diff$diffgenes))
      }
    }
    }else{
      cg = m$symbol[m$symbol %in% my_genes]
    }
  heatmap = draw_heatmap(exp[cg,],
                         group_list,
                         ...)
  return(heatmap)
}


##' multi_deg_all
##'
##' do diffiencial analysis according to exprission set and group information
##'
##' @inheritParams draw_pca
##' @inheritParams draw_volcano
##' @inheritParams draw_heatmap
##' @inheritParams draw_heatmap2
##' @inheritParams multi_deg
##' @param color_volcano color for volcano
##' @return a list with deg data.frame, volcano plot and a list with DEGs.
##' @author Xiaojie Sun
##' @importFrom patchwork wrap_plots
##' @importFrom stringr str_split
##' @importFrom dplyr union_all
##' @importFrom patchwork plot_layout
##' @export
##' @examples
##' \donttest{
##' gse = "GSE474"
##' geo = geo_download(gse,destdir=tempdir(),by_annopbrobe = FALSE)
##' geo$exp[1:4,1:4]
##' geo$exp=log2(geo$exp+1)
##' group_list=ifelse(stringr::str_detect(geo$pd$title,"MObese"),"MObese",
##' ifelse(stringr::str_detect(geo$pd$title,"NonObese"),"NonObese","Obese"))
##' group_list=factor(group_list,levels = c("NonObese","Obese","MObese"))
##' find_anno(geo$gpl)
##' ids = AnnoProbe::idmap(geo$gpl,destdir = tempdir())
##' dcp = multi_deg_all(geo$exp,
##' group_list,ids,adjust = FALSE)
##' dcp[[3]]
##' }
##' @seealso
##' \code{\link{geo_download}};\code{\link{draw_volcano}};\code{\link{draw_venn}}

multi_deg_all <- function(exp,
                          group_list,
                          ids,
                          pkg = 4,
                          symmetry = TRUE,
                          my_genes = NULL,
                          show_rownames = FALSE,
                          cluster_cols = TRUE,
                          ...
                          ) {
  deg = multi_deg(
    exp,
    group_list,
    ids,
    ...
  )
  #exp = data.frame(exp)
  #exp = exp[match(deg[[1]]$probe_id,rownames(exp)),]
  cgs = get_cgs(deg)
  volcano_plot = draw_volcano2(deg,
                               pkg = pkg,
                               symmetry = symmetry)
  pca_plot = draw_pca(exp,group_list)
  heatmap = draw_heatmap2(exp,group_list,deg,my_genes,show_rownames = show_rownames,cluster_cols = cluster_cols)
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
  result = list(
    deg = deg,
    cgs = cgs,
    plots = wrap_plots(plotlist) +
      plot_layout(design = layout) +
      plot_layout(guides = 'collect')
  )
  diffprobes = lapply(cgs,function(x)x$diff$diffprobes)
  message(paste0(length(union_all(diffprobes))," DEGs in all,",length(intersect_all(diffprobes))," DEGs in common."))
  return(result)
  }
