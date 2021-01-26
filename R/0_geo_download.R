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
  if((by_annopbrobe) & !require(AnnoProbe)) stop("you must install AnnoProbe first by devtools::install_github('jmzeng1314/AnnoProbe')")
  if(by_annopbrobe){
    if((!file.exists(paste0(gse,"_eSet.Rdata")))) geoChina(gse)
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
