##' multi_deg
##'
##' do diffiential analysis according to expression set and group information
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
##' \dontrun{
##' gse = "GSE474"
##' geo = geo_download(gse,destdir=tempdir())
##' geo$exp[1:4,1:4]
##' geo$exp=log2(geo$exp+1)
##' group_list=ifelse(stringr::str_detect(geo$pd$title,"MObese"),
##' "MObese",ifelse(stringr::str_detect(geo$pd$title,"NonObese"),
##' "NonObese","Obese"))
##' group_list=factor(group_list,levels = c("NonObese","Obese","MObese"))
##' find_anno(geo$gpl)
##' ids <- AnnoProbe::idmap(geo$gpl,destdir = tempdir())
##' deg = multi_deg(geo$exp,group_list,ids,adjust = FALSE)
##' names(deg)
##' head(deg[[1]])
##' head(deg[[2]])
##' head(deg[[3]])
##' }
##' @seealso
##' \code{\link{get_deg}};\code{\link{multi_deg_all}}

multi_deg <- function(exp,
                      group_list,
                      ids,
                      logFC_cutoff = 1,
                      pvalue_cutoff = 0.05,
                      adjust = FALSE,
                      species = "human",
                      entriz = TRUE) {
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
    design=stats::model.matrix(~0+group_list)
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
    design=stats::model.matrix(~0+group_list)
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
    design=stats::model.matrix(~0+group_list)
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
    deg[[i]] <- mutate(deg[[i]],probe_id=rownames(deg[[i]]))
    ids = stats::na.omit(ids)
    if(is.character(ids$probe_id)) ids$probe_id = as.character(ids$probe_id)
    ids = ids[!duplicated(ids$symbol),]
    ids = ids[!duplicated(ids$probe_id),]
    deg[[i]] <- inner_join(deg[[i]],ids,by="probe_id")
    if(adjust){
      k1 = (deg[[i]]$adj.P.Val < pvalue_cutoff)&(deg[[i]]$logFC < -logFC_cutoff)
      k2 = (deg[[i]]$adj.P.Val < pvalue_cutoff)&(deg[[i]]$logFC > logFC_cutoff)
    }else{
      k1 = (deg[[i]]$P.Value < pvalue_cutoff)&(deg[[i]]$logFC < -logFC_cutoff)
      k2 = (deg[[i]]$P.Value < pvalue_cutoff)&(deg[[i]]$logFC > logFC_cutoff)
    }
    change = ifelse(k1,"down",ifelse(k2,"up","stable"))
    deg[[i]] <- mutate(deg[[i]],change)
    if(entriz){
      if(species == "human"){
        if(!requireNamespace("org.Hs.eg.db",quietly = TRUE)) {
          stop("Package \"org.Hs.eg.db\" needed for this function to work.
         Please install it by BiocManger::install('org.Hs.eg.db')",call. = FALSE)
        }
        or = org.Hs.eg.db::org.Hs.eg.db
      }
      if(species == "mouse"){
        if(!requireNamespace("org.Mm.eg.db",quietly = TRUE)) {
          stop("Package \"org.Mm.eg.db\" needed for this function to work.
         Please install it by BiocManger::install('org.Mm.eg.db')",call. = FALSE)
        }
        or = org.Mm.eg.db::org.Mm.eg.db
      }
      if(species == "rat"){
        if(!requireNamespace("org.Rn.eg.db",quietly = TRUE)) {
          stop("Package \"org.Rn.eg.db\" needed for this function to work.
         Please install it by BiocManger::install('org.Rn.eg.db')",call. = FALSE)
        }
        or = org.Rn.eg.db::org.Rn.eg.db
      }
      s2e <- bitr(unique(deg[[i]]$symbol), fromType = "SYMBOL",
                  toType = c( "ENTREZID"),
                  OrgDb = or)
      s2e <- s2e[!duplicated(s2e$SYMBOL),]
      deg[[i]] <- inner_join(deg[[i]],s2e,by=c("symbol"="SYMBOL"))
    }
  }
  return(deg)
}
