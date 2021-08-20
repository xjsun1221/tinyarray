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
##' group_list=ifelse(str_detect(geo$pd$title,"MObese"),"MObese",
##' ifelse(str_detect(geo$pd$title,"NonObese"),"NonObese","Obese"))
##' group_list=factor(group_list,levels = c("NonObese","Obese","MObese"))
##' find_anno(geo$gpl)
##' library(hgu133a.db)
##' ids <- toTable(hgu133aSYMBOL)
##' deg = multi_deg(geo$exp,group_list,ids,adjust = FALSE)
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
##' group_list=ifelse(str_detect(geo$pd$title,"MObese"),"MObese",
##' ifelse(str_detect(geo$pd$title,"NonObese"),"NonObese","Obese"))
##' group_list=factor(group_list,levels = c("NonObese","Obese","MObese"))
##' find_anno(geo$gpl)
##' library(hgu133a.db)
##' ids <- toTable(hgu133aSYMBOL)
##' deg = multi_deg(geo$exp,group_list,ids,adjust = FALSE)
##' draw_volcano2(deg)
##' draw_volcano2(deg,color = c("darkgreen","grey","darkred"))
##' @seealso
##' \code{\link{geo_download}};\code{\link{draw_volcano}};\code{\link{draw_venn}}

draw_volcano2 = function(deg,
                         pkg=4,
                         lab,
                         pvalue_cutoff=0.05,
                         logFC_cutoff=1,
                         adjust=F,
                         symmetry=T,
                         color = c("blue", "grey", "red")
){
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
                              color = color,
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
##' @param heat_union logical ,use union or intersect DEGs for heatmap
##' @param heat_id id of heatmap,1 for all DEGs,2 for head and tail,3 for top n DEGs
##' @param gene_number how many DEGs will heatmap show .
##' @return a heatmap plot according to \code{exp} and grouped by \code{group}.
##' @author Xiaojie Sun
##' @importFrom pheatmap pheatmap
##' @export
##' @examples
##' gse = "GSE474"
##' geo = geo_download(gse)
##' geo$exp[1:4,1:4]
##' geo$exp=log2(geo$exp+1)
##' library(stringr)
##' group_list=ifelse(str_detect(geo$pd$title,"MObese"),"MObese",
##' ifelse(str_detect(geo$pd$title,"NonObese"),"NonObese","Obese"))
##' group_list=factor(group_list,levels = c("NonObese","Obese","MObese"))
##' find_anno(geo$gpl)
##' library(hgu133a.db)
##' ids <- toTable(hgu133aSYMBOL)
##' deg = multi_deg(geo$exp,group_list,ids,adjust = FALSE)
##' draw_heatmap2(geo$exp,group_list,deg)
##' @seealso
##' \code{\link{draw_pca}};\code{\link{draw_volcano}};\code{\link{draw_venn}}

draw_heatmap2 <- function(exp,
                          group_list,
                          deg,
                          heat_union = TRUE,
                          heat_id=1,
                          gene_number=200,
                          show_rownames = FALSE,
                          scale_before = FALSE,
                          n_cutoff = 3,
                          cluster_cols = TRUE,
                          annotation_legend=F,
                          legend = FALSE,
                          color = grDevices::colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
                          color_an = c("#92C5DE", "#F4A582", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3",
                                       "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3")
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
  np2 = np[order(apply(np, 1, stats::mad),decreasing = T),]
  n = switch(heat_id,
             all = np,
             ht  = rbind(utils::head(np2,gene_number),
                         utils::tail(np2,gene_number)),
             top = utils::head(np,gene_number))
  heatmap = draw_heatmap(n,
                         group_list,
                         legend = legend,
                         show_rownames = show_rownames,
                         scale_before = scale_before,
                         n_cutoff = n_cutoff,
                         cluster_cols = cluster_cols,
                         annotation_legend = annotation_legend,
                         color = color,
                         color_an = color_an
  )
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
##' gse = "GSE474"
##' geo = geo_download(gse)
##' geo$exp[1:4,1:4]
##' geo$exp=log2(geo$exp+1)
##' library(stringr)
##' group_list=ifelse(str_detect(geo$pd$title,"MObese"),"MObese",
##' ifelse(str_detect(geo$pd$title,"NonObese"),"NonObese","Obese"))
##' group_list=factor(group_list,levels = c("NonObese","Obese","MObese"))
##' find_anno(geo$gpl)
##' library(hgu133a.db)
##' ids <- toTable(hgu133aSYMBOL)
##' dcp = multi_deg_all(geo$exp,
##' group_list,ids,adjust = FALSE)
##' dcp[[3]]
##' @seealso
##' \code{\link{geo_download}};\code{\link{draw_volcano}};\code{\link{draw_venn}}

multi_deg_all <- function(exp,
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
                          color_volcano = c("blue","grey","red")) {
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
  volcano_plot = draw_volcano2(deg,
                               pkg = pkg,
                               pvalue_cutoff = pvalue_cutoff,
                               logFC_cutoff = logFC_cutoff,
                               adjust = adjust,
                               symmetry = symmetry,
                               color = color_volcano)
  pca_plot = draw_pca(exp,group_list)
  heatmap = draw_heatmap2(exp,group_list,
                          deg,
                          show_rownames = show_rownames,
                          heat_id=heat_id,
                          gene_number=gene_number,
                          scale_before = scale_before,
                          n_cutoff = n_cutoff,
                          cluster_cols = cluster_cols,
                          annotation_legend=annotation_legend
                          )
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
  print(paste0(length(union_all(diffprobes))," DEGs in all,",length(intersect_all(diffprobes))," DEGs in common."))
  return(result)
  }
