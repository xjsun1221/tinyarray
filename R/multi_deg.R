##' multi_deg
##'
##' do diffiencial analysis according to exprission set and group information
##'
##' @inheritParams get_deg
##' @param ids a data.frame with 2 columns,including probe_id and symbol
##' @return a deg data.frame with 10 columns
##' @author Xiaojie Sun
##' @importFrom limma lmFit
##' @importFrom limma eBayes
##' @importFrom limma topTable
##' @importFrom limma makeContrasts
##' @importFrom limma contrasts.fit
##' @importFrom clusterProfiler bitr
##' @importFrom dplyr mutate
##' @importFrom dplyr inner_join
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
##' @seealso
##' \code{\link{geo_download}};\code{\link{draw_volcano}};\code{\link{draw_venn}}

multi_deg <- function(exp,
                      group_list,
                      ids,
                      logFC_cutoff = 1,
                      pvalue_cutoff = 0.05,
                      adjust = F,
                      entriz = T) {
  p1 <-  all(apply(exp,2,is.numeric))
  if(!p1) stop("exp must be a numeric matrix")
  p2  <-  (sum(!duplicated(group_list)) > 1)
  if(!p2) stop("group_list must more than 1")
  p3 <- is.factor(group_list)
  if(!p3) stop("group_list must be a factor")
  if(ncol(exp)!=length(group_list))stop("wrong group_list or exp")
  if(ncol(ids)!=2)stop("wrong ids pramater,it should be a data.frame with probe_id and symbol")
  colnames(ids) = c("probe_id","symbol")
  px <-  levels(group_list)
  if(length(px)==3){
    design=model.matrix(~0+group_list)
    colnames(design)=c("x1","x2","x3")
    px <-  levels(group_list)
    contrast.matrix <- makeContrasts("x2-x1",
                                     "x3-x1",
                                     "x3-x2",
                                     levels=design)
    rownames(contrast.matrix) = levels(group_list)
    colnames(contrast.matrix) = c(paste0(px[2],"-",px[1]),
                                  paste0(px[3],"-",px[1]),
                                  paste0(px[3],"-",px[2]))
    colnames(design) = levels(group_list)
  }else if(length(px)==4){
    design=model.matrix(~0+group_list)
    colnames(design)=c("x1","x2","x3","x4")
    px <-  levels(group_list)
    contrast.matrix <- makeContrasts("x2-x1",
                                     "x3-x1",
                                     "x4-x1",
                                     levels=design)
    rownames(contrast.matrix) = levels(group_list)
    colnames(contrast.matrix) = c(paste0(px[2],"-",px[1]),
                                  paste0(px[3],"-",px[1]),
                                  paste0(px[4],"-",px[1]))
    colnames(design) = levels(group_list)
  }else if(length(px)==5){
    design=model.matrix(~0+group_list)
    colnames(design)=c("x1","x2","x3","x4","x5")
    px <-  levels(group_list)
    contrast.matrix <- makeContrasts("x2-x1",
                                     "x3-x1",
                                     "x4-x1",
                                     "x5-x1",
                                     levels=design)
    rownames(contrast.matrix) = levels(group_list)
    colnames(contrast.matrix) = c(paste0(px[2],"-",px[1]),
                                  paste0(px[3],"-",px[1]),
                                  paste0(px[4],"-",px[1]),
                                  paste0(px[5],"-",px[1]))
    colnames(design) = levels(group_list)
  }
  fit <- lmFit(exp, design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  if(length(px)>3) n = 1:(length(px)-1) else{n = 1:length(px)}
  deg = lapply(n, function(i){
    topTable(fit2,coef = i,number = Inf)
  })
  names(deg) = colnames(contrast.matrix)

  for(i in 1:length(deg)){
    #1.加probe_id列
    deg[[i]] <- mutate(deg[[i]],probe_id=rownames(deg[[i]]))
    #2.id转换
    ids = na.omit(ids[match(names(sort(apply(exp, 1, mad),decreasing = T)),
                             ids$probe_id),])
    ids = ids[!duplicated(ids$symbol),]
    ids = ids[!duplicated(ids$probe_id),]
    deg[[i]] <- inner_join(deg[[i]],ids,by="probe_id")
    #deg[[i]] <- deg[[i]][!duplicated(deg[[i]]$symbol),]
    #3.change
    if(adjust){
      k1 = (deg[[i]]$adj.P.Val < pvalue_cutoff)&(deg[[i]]$logFC < -logFC_cutoff)
      k2 = (deg[[i]]$adj.P.Val < pvalue_cutoff)&(deg[[i]]$logFC > logFC_cutoff)
    }else{
      k1 = (deg[[i]]$P.Value < pvalue_cutoff)&(deg[[i]]$logFC < -logFC_cutoff)
      k2 = (deg[[i]]$P.Value < pvalue_cutoff)&(deg[[i]]$logFC > logFC_cutoff)
    }
    change = ifelse(k1,"down",ifelse(k2,"up","stable"))
    deg[[i]] <- mutate(deg[[i]],change)
    #4.加ENTREZID列，后面富集分析要用
    if(entriz){
    s2e <- bitr(unique(deg[[i]]$symbol), fromType = "SYMBOL",
                toType = c( "ENTREZID"),
                OrgDb = org.Hs.eg.db)
    s2e <- s2e[!duplicated(s2e$SYMBOL),]
    deg[[i]] <- inner_join(deg[[i]],s2e,by=c("symbol"="SYMBOL"))
    }
  }
  return(deg)

}
