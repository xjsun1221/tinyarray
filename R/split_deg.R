##' get_cgs
##'
##' extract DEGs from deg data.frame
##'
##' @inheritParams draw_volcano
##' @return a list with upgenes,downgenes,diffgenes.
##' @author Xiaojie Sun
##' @export
##' @examples
##' #two group
##' gse = "GSE42872"
##' geo = geo_download(gse)
##' group_list = rep(c("A","B"),each = 3)
##' ids = AnnoProbe::idmap('GPL6244')
##' deg = get_deg(geo$exp,group_list,ids)
##' cgs = get_cgs(deg)
##' #mutigroup
##' gse = "GSE474"
##' geo = geo_download(gse)
##' geo$exp[1:4,1:4]
##' geo$exp=log2(geo$exp+1)
##' library(stringr)
##' group_list=ifelse(str_detect(geo$pd$title,"MObese"),"MObese",ifelse(str_detect(geo$pd$title,"NonObese"),"NonObese","Obese"))
##' group_list=factor(group_list,levels = c("NonObese","Obese","MObese"))
##' find_anno(geo$gpl)
##' ids <- toTable(hgu133aSYMBOL)
##' deg = multi_deg(geo$exp,group_list,ids,adjust = F)
##' cgs = get_cgs(deg)
##' @seealso
##' \code{\link{geo_download}};\code{\link{draw_volcano}};\code{\link{draw_venn}}


get_cgs <- function(deg){
  if(!is.list(deg) & is.data.frame(deg))stop("deg is a data.frame or list returned by limma")
  cgs = list()
  if(is.data.frame(deg)) deg = list(deg = deg)
  for(i in 1:length(deg)){
    cgs[[i]] = list(up = data.frame(upgenes =deg[[i]]$symbol[deg[[i]]$change=="up"],
                                         upprobes = deg[[i]]$probe_id[deg[[i]]$change=="up"],
                                    stringsAsFactors = F),
                    down = data.frame(downgenes = deg[[i]]$symbol[deg[[i]]$change=="down"],
                                           downprobes = deg[[i]]$probe_id[deg[[i]]$change=="down"],
                                      stringsAsFactors = F),
                    diff = data.frame(diffgenes = deg[[i]]$symbol[deg[[i]]$change!="stable"],
                                           diffprobes = deg[[i]]$probe_id[deg[[i]]$change!="stable"],
                                      stringsAsFactors = F)
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
##' #two group
##' gse = "GSE42872"
##' geo = geo_download(gse)
##' group_list = rep(c("A","B"),each = 3)
##' ids = AnnoProbe::idmap('GPL6244')
##' deg = get_deg(geo$exp,group_list,ids)
##' draw_volcano2(deg)
##' #multigroup
##' gse = "GSE474"
##' geo = geo_download(gse)
##' geo$exp[1:4,1:4]
##' geo$exp=log2(geo$exp+1)
##' library(stringr)
##' group_list=ifelse(str_detect(geo$pd$title,"MObese"),"MObese",ifelse(str_detect(geo$pd$title,"NonObese"),"NonObese","Obese"))
##' group_list=factor(group_list,levels = c("NonObese","Obese","MObese"))
##' find_anno(geo$gpl)
##' ids <- toTable(hgu133aSYMBOL)
##' deg = multi_deg(geo$exp,group_list,ids,adjust = F)
##' draw_volcano2(deg)
##' @seealso
##' \code{\link{geo_download}};\code{\link{draw_volcano}};\code{\link{draw_venn}}

draw_volcano2 = function(deg,
                         pkg=4,
                         lab,
                         pvalue_cutoff=0.05,
                         logFC_cutoff=1,
                         adjust=F,
                         symmetry=T){
  if(!is.list(deg) & is.data.frame(deg))stop("deg is a data.frame or list returned by limma")
  if(is.data.frame(deg)) deg = list(deg = deg)
  volcano_plots <- lapply(1:length(deg),
                          function(k) {
                            draw_volcano(
                              deg[[k]] ,
                              pkg = pkg,
                              lab = names(deg)[k],
                              pvalue_cutoff = pvalue_cutoff,
                              logFC_cutoff = logFC_cutoff,
                              adjust = adjust,
                              symmetry = symmetry
                            )
                          })
  if(is.data.frame(deg)) {
    volcano_plots = volcano_plots[[1]]
  }else{
    volcano_plots = wrap_plots(volcano_plots)+
      plot_layout(design = paste(rep(LETTERS[1:length(deg)]),collapse = "")) +
      plot_layout(guides = 'collect')
    }
  return(wrap_plots(volcano_plots))
}

##' draw heatmap plots
##'
##' print heatmap plots for expression matrix and group by group_list paramter
##'
##' @inheritParams draw_volcano
##' @inheritParams draw_heatmap
##' @inheritParams draw_pca
##' @param heat_union logical , if TRUE ,use union or intersect DEGs for heatmap
##' @param heat_id id of heatmap,1 for all DEGs,2 for head and tail,3 for top n DEGs
##' @param gene_number how many DEGs will heatmap show .
##' @return a heatmap plot according to \code{exp} and grouped by \code{group}.
##' @author Xiaojie Sun
##' @importFrom pheatmap pheatmap
##' @importFrom ggplotify as.ggplot
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
##' deg = multi_deg(geo$exp,group_list,ids,adjust = F)
##' draw_heatmap2(geo$exp,group_list,deg)
##' @seealso
##' \code{\link{draw_pca}};\code{\link{draw_volcano}};\code{\link{draw_venn}}


draw_heatmap2 <- function(exp,
                          group_list,
                          deg,
                          heat_union = T,
                          heat_id=1,
                          gene_number=200,
                          scale_before = F,
                          n_cutoff = 3,
                          cluster_cols = T,
                          annotation_legend=F,
                          color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)
                          ){
  cgs = get_cgs(deg)
  if(length(cgs)==1){
    cg = cgs[[1]]$diff$diffprobes
    cgup = cgs[[1]]$up$upprobes
    cgdown = cgs[[1]]$down$downprobes
  }else{
    if(heat_union){
      cg = union_all(lapply(cgs,function(x)x$diff$diffprobes))
      cgup = union_all(lapply(cgs,function(x)x$up$upprobes))
      cgdown = union_all(lapply(cgs,function(x)x$down$downprobes))
    }else{
      cg = intersect_all(lapply(cgs,function(x)x$diff$diffprobes))
      cgup = intersect_all(lapply(cgs,function(x)x$up$upprobes))
      cgdown = intersect_all(lapply(cgs,function(x)x$down$downprobes))
    }
  }
  cd = gene_number > length(cgup) | gene_number > length(cgdown)
  er = paste0("gene_number must less than ",min(c(length(cgup),length(cgdown))))
  if(heat_id==2 & cd) stop(er)
  if(heat_id==3 & gene_number> length(cg))stop(paste0("gene_number must less than ",length(cg)))
  np = exp[cg,]
  np2 = np[order(apply(np, 1, mad),decreasing = T),]
  n = switch(heat_id,
             all = np,
             ht  = rbind(head(np2,gene_number),
                         tail(np2,gene_number)),
             top = head(np,gene_number))
  heatmap = draw_heatmap(n,
                         group_list,
                         scale_before = scale_before,
                         n_cutoff = n_cutoff,
                         cluster_cols = cluster_cols,
                         annotation_legend = annotation_legend,
                         color = color
  )
}
