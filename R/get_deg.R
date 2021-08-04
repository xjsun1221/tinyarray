##' get_deg
##'
##' do diffiencial analysis according to exprission set and group information
##'
##' @inheritParams draw_pca
##' @inheritParams draw_volcano
##' @param ids a data.frame with 2 columns,including probe_id and symbol
##' @return a deg data.frame with 10 columns
##' @author Xiaojie Sun
##' @importFrom limma lmFit
##' @importFrom limma eBayes
##' @importFrom limma topTable
##' @importFrom clusterProfiler bitr
##' @importFrom dplyr mutate
##' @importFrom dplyr inner_join
##' @export
##' @examples
##' gse = "GSE42872"
##' geo = geo_download(gse)
##' Group = rep(c("control","treat"),each = 3)
##' Group = factor(Group)
##' find_anno(geo$gpl,install = T)
##' ids <- toTable(hugene10sttranscriptclusterSYMBOL)
##' deg = get_deg(geo$exp,Group,ids)
##' head(deg)
##' @seealso
##' \code{\link{geo_download}};\code{\link{draw_volcano}};\code{\link{draw_venn}}

get_deg <- function(exp,
                    group_list,
                    ids,
                    logFC_cutoff=1,
                    pvalue_cutoff=0.05,
                    adjust = F,
                    entriz = T) {
  p3 <- is.factor(group_list)
  if(!p3) {
    group_list = factor(group_list)
    warning("group_list was covert to factor")
  }
  if(ncol(exp)!=length(group_list))stop("wrong group_list or exp")
  if(ncol(ids)!=2)stop("wrong ids pramater,it should be a data.frame with probe_id and symbol")
  colnames(ids) = c("probe_id","symbol")

  design=model.matrix(~group_list)
  fit=lmFit(exp,design)
  fit=eBayes(fit)
  deg=topTable(fit,coef=2,number = Inf)

  #为deg数据框添加几列
  #1.加probe_id列，把行名变成一列
  if("ID" %in% colnames(deg)){
  deg <- mutate(deg,probe_id=deg$ID)
}else{
  deg <- mutate(deg,probe_id=rownames(deg))
}
  head(deg)
  #2.加symbol列，火山图要用
  ids = ids[!duplicated(ids$symbol),]
  deg <- inner_join(deg,ids,by="probe_id")
  head(deg)
  #按照symbol列去重复
  #deg <- deg[!duplicated(deg$symbol),]
  #3.加change列,标记上下调基因
  if(adjust){
    k1 = (deg$adj.P.Val < pvalue_cutoff)&(deg$logFC < -logFC_cutoff)
    k2 = (deg$adj.P.Val < pvalue_cutoff)&(deg$logFC > logFC_cutoff)
  }else{
    k1 = (deg$P.Value < pvalue_cutoff)&(deg$logFC < -logFC_cutoff)
    k2 = (deg$P.Value < pvalue_cutoff)&(deg$logFC > logFC_cutoff)
  }

  change = ifelse(k1,
                  "down",
                  ifelse(k2,
                         "up",
                         "stable"))
  deg <- mutate(deg,change)
  #4.加ENTREZID列，用于富集分析（symbol转entrezid，然后inner_join）

  if(entriz){
    s2e <- bitr(deg$symbol,
                fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = org.Hs.eg.db::org.Hs.eg.db)#人类
    #其他物种http://bioconductor.org/packages/release/BiocViews.html#___OrgDb
    deg <- inner_join(deg,s2e,by=c("symbol"="SYMBOL"))
    deg <- deg[!duplicated(deg$symbol),]
  }
  return(deg)
}
