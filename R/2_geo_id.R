##' geo_download
##'
##' download gse data and get informations
##'
##' @param gse gse assession number
##' @param by_annopbrobe download data by geoquery or annoprobe
##' @param simpd get simplified pdata,drop out columns with all same values
##' @param colon_remove whether to remove duplicated columns with colons
##' @param destdir	 The destination directory for data downloads.
##' @param n For data with more than one ExpressionSet, specify which one to analyze
##' @return a list with exp,pd and gpl
##' @author Xiaojie Sun
##' @export
##' @examples
##' \dontrun{
##' gse = "GSE42872"
##' a = geo_download(gse,destdir=tempdir())
##' }
##' @seealso
##' \code{\link{find_anno}}

geo_download <-  function(gse,by_annopbrobe = TRUE,
                          simpd = TRUE,colon_remove = FALSE,
                          destdir = getwd(),n = 1){

  if(!requireNamespace("Biobase",quietly = TRUE)) {
    stop("Package \"Biobase\" needed for this function to work.
         Please install it by BiocManager::install('Biobase')",call. = FALSE)
  }
  if((!by_annopbrobe) & !requireNamespace("GEOquery",quietly = TRUE)) {
    stop("Package \"GEOquery\" needed for this function to work.
         Please install it by BiocManager::install('GEOquery')",call. = FALSE)
  }
  if((by_annopbrobe) & !requireNamespace("AnnoProbe",quietly = TRUE)) {
    stop("Package \"Biobase\" needed for this function to work.
         Please install it by install.packages('AnnoProbe')",call. = FALSE)
  }
  if(by_annopbrobe){
    if(!file.exists(paste0(destdir,"/",gse,"_eSet.Rdata"))){
      eSet <- tryCatch({AnnoProbe::geoChina(gse, destdir = destdir)
      },error = function(e){555})

      if(eSet == 555){
        warning("This data is not indexed by AnnoProbe, downloaded by GEOquery")
        eSet <- GEOquery::getGEO(gse,destdir = destdir,getGPL = FALSE)
        gset = eSet
        save(gset,file = paste0(destdir,"/",gse,"_eSet.Rdata"))
      }
    }else{
      suppressWarnings(load(paste0(destdir,"/",gse,"_eSet.Rdata")))
      eSet = gset
      rm(gset)
    }
  }else{
    eSet <- GEOquery::getGEO(gse,destdir = destdir,getGPL = FALSE)
  }
  if(length(n)!=1) stop("only one ExpresssionSet can be analyzed")
  if(length(eSet)==1 & n!=1) {
    n = 1
    warning("this data only have one ExpresssionSet object")
  }

  exp <- Biobase::exprs(eSet[[n]])
  pd <- Biobase::pData(eSet[[n]])
  if(simpd){
    colname <- vector("character")
    count <- vector("integer")
    for (i in 1:ncol(pd)) {
      colname[i] = colnames(pd)[[i]]
      count[i] = nrow(pd[!duplicated(pd[, i]), ])
    }
    df <- data.frame(colname, count) |>
      dplyr::arrange(desc(count)) |> dplyr::filter(count >1)
    pd <-  pd[,df$colname]
    pd <- pd[,!(colnames(pd) %in% c("geo_accession", "supplementary_file"))]
    if(colon_remove){
      pd = pd[,!apply(pd, 2, function(x){all(stringr::str_detect(x,": |https|www")|is.na(x)|x=="")})]
      colnames(pd) = stringr::str_remove(colnames(pd),":ch1")
    }
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
  gpl <- eSet[[n]]@annotation
  pd2 = apply(pd,2,as.character) |> as.data.frame()
  rownames(pd2) = rownames(pd)
  re = list(exp=exp,pd=pd2,gpl=gpl)
  if(is.null(dim(exp)) | nrow(exp)==0){
    warning("exp is empty")
  } else if (any(is.na(exp)|is.nan(exp))) {
    warning("NA or NAN values detected")
  }else if (any(exp<0)) {
    warning("nagtive values detected")
  } else{
    message(paste(nrow(exp),"probes,",
                  ncol(exp),"samples",
                  "from",min(exp),
                  "to",max(exp)))}
  return(re)
}

##' find annotation package or files
##'
##' find gpl annotation package or files
##'
##' @param gpl a gpl accession
##' @param install whether to install and library the package
##' @param update whether to update the package
##' @return a list with deg data.frame, volcano plot and a list with DEGs.
##' @author Xiaojie Sun
##' @importFrom stringr str_remove_all
##' @importFrom stringr str_to_upper
##' @importFrom BiocManager install
##' @export
##' @examples
##' find_anno("GPL570")
##' @seealso
##' \code{\link{geo_download}}

find_anno <-function(gpl,install = FALSE,update = FALSE){
  gpl = str_to_upper(gpl)
  if(!any(pkg_all$gpl==gpl)) {
    if(gpl %in% setdiff(exists_anno_list,pkg_all$gpl)){
      ml1 = str_remove_all(paste0("`ids <- AnnoProbe::idmap\\(","\\'",gpl,"\\'","\\)`"),"\\\\")
      message(paste0("no annotation packages avliable,please use ",ml1))
      message("if you get error by this code ,please try different `type` parameters")
    }else{
      message("no annotation avliable in Bioconductor and AnnoProbe")
    }
  }else {
    qz = pkg_all$bioc_package[pkg_all$gpl== gpl]
    ml1 = str_remove_all(paste0("`ids <- AnnoProbe::idmap\\(","\\'",gpl,"\\'","\\)`"),"\\\\")
    ml2 = str_remove_all(paste0("`library\\(",qz,".db","\\)",";","ids <- toTable\\(",qz,"SYMBOL\\)`"),"\\\\")
    if(install){
      if(!suppressMessages(requireNamespace(paste0(qz,".db")))){
        BiocManager::install(paste0(qz,".db"),update = update,ask = FALSE)
        suppressMessages(requireNamespace(paste0(qz,".db")))
      }
    }
    if(!(gpl %in% exists_anno_list)) {
      message(paste0(ml2," is avaliable"))
    }else {
      message(paste0(ml2," and ",ml1 ," are both avaliable"))
      message("if you get error by idmap, please try different `type` parameters")
    }
  }
}
utils::globalVariables(c("pkg_all","exists_anno_list","gset"))
