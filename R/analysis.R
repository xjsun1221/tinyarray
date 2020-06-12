##' geo_download
##'
##' download gse data and get informations
##'
##' @param gse gse assession number
##' @param by_annopbrobe getGEO or geoChina
##' @param simpd get simplified pdata,drop out columns with all same values
##' @return a list with exp,pd and gpl
##' @author Xiaojie Sun
##' @importFrom GEOquery getGEO
##' @importFrom Biobase exprs
##' @importFrom Biobase pData
##' @importFrom AnnoProbe geoChina
##' @importFrom dplyr arrange
##' @importFrom dplyr filter
##' @importFrom dplyr %>%
##' @export
##' @examples
##' gse = "GSE42872"
##' geo_download(gse)
##' geo_download(gse,by_annopbrobe = F)
##' @seealso
##' \code{\link{simpd}};\code{\link{draw_volcano}};\code{\link{draw_venn}}

geo_download <-  function(gse,by_annopbrobe = T,simpd=T){
  if((by_annopbrobe = T) & !require(AnnoProbe)) stop("you must install AnnoProbe first by devtools::install_github('jmzeng1314/AnnoProbe')")
  if(by_annopbrobe){
    if(!file.exists(paste0(gse,"_eSet.Rdata"))) geoChina(gse)
    load(paste0(gse,"_eSet.Rdata"))
    eSet <- gset
    rm(gset)
  }else{
    eSet <- getGEO(gse,
                   destdir = '.',
                   getGPL = F)
  }
  #(1)提取表达矩阵exp
  exp <- exprs(eSet[[1]])
  #(2)提取临床信息
  pd <- pData(eSet[[1]])
  if(simpd){
    colname <- vector("character")
    count <- vector("integer")
    for (i in 1:ncol(pd)) {
      colname[i] = colnames(pd)[[i]]
      count[i] = nrow(pd[!duplicated(pd[, i]), ])
    }
    df <- data.frame(colname, count,stringsAsFactors = F) %>% arrange(desc(count)) %>% dplyr::filter(count >1)
    pd <-  pd[,df$colname]
  }
  p1 = identical(rownames(pd),colnames(exp))
  p2 = all(rownames(pd) %in% colnames(exp) & colnames(exp) %in% rownames(pd))
  if(!p1) {
    exp = exp[,match(rownames(pd),colnames(exp))]
    if(!p2) {
      exp = exp[,intersect(rownames(pd),colnames(exp))]
      pd = pd[intersect(rownames(pd),colnames(exp)),]
    }
  }
  #(3)提取芯片平台编号
  gpl <- eSet[[1]]@annotation
  re = list(exp=exp,pd=pd,gpl=gpl)
  return(re)
}

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
##' geo_download(gse)
##' geo_download(gse,by_annopbrobe = F)
##' @seealso
##' \code{\link{geo_download}};\code{\link{draw_volcano}};\code{\link{draw_venn}}

get_deg <- function(exp,
                    group_list,
                    ids,
                    logFC_cutoff=1,
                    pvalue_cutoff=0.05,
                    adjust = F,
                    entriz = T) {
  if(ncol(exp)!=length(group_list))stop("wrong group_list or exp")
  if(ncol(ids)!=2)stop("wrong ids pramater,it should be a data.frame with probe_id and symbol")
  colnames(ids) = c("probe_id","symbol")

  design=model.matrix(~group_list)
  fit=lmFit(exp,design)
  fit=eBayes(fit)
  deg=topTable(fit,coef=2,number = Inf)

  #为deg数据框添加几列
  #1.加probe_id列，把行名变成一列
  deg <- mutate(deg,probe_id=rownames(deg))
  head(deg)
  #2.加symbol列，火山图要用
  deg <- inner_join(deg,ids,by="probe_id")
  head(deg)
  #按照symbol列去重复
  deg <- deg[!duplicated(deg$symbol),]
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
  }
  return(deg)
}
