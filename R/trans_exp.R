##' trans_exp
##'
##' transform rownames of tcga or tcga_gtex expression set from gdc or xena,from ensamble id to gene symbol
##'
##' @param exp tcga or tcga_gtex expression set from gdc or xena
##' @param mrna_only only keep mrna rows in result
##' @param lncrna_only only keep lncrna rows in result
##' @param gtex logical,whether including gtex data
##' @return a transformed expression set with symbol
##' @author Xiaojie Sun
##' @importFrom stringr str_detect
##' @importFrom stringr str_remove
##' @importFrom dplyr inner_join
##' @export
##' @examples
##' exp = matrix(rnorm(1000),ncol = 10)
##' rownames(exp) = sample(mRNA_annov23$gene_id,100)
##' colnames(exp) = c(paste0("TCGA",1:5),paste0("GTEX",1:5))
##' k  = trans_exp(exp)

##' @seealso
##' \code{\link{simpd}};\code{\link{draw_volcano}};\code{\link{draw_venn}}

trans_exp = function(exp,mrna_only = F,lncrna_only = F,gtex = F){
  k00 = any(str_detect(colnames(exp),"TCGA"))
  if(!k00)warning("this expression set probably not from TCGA,please ensure it")
  k0 = any(str_detect(colnames(exp),"GTEX"))
  kd = any(str_detect(rownames(exp),"\\."))
  if((!(k0|gtex))){
    lanno = lnc_anno
    manno = mRNA_anno
  }else if(k00){
    lanno = lnc_annov23
    manno = mRNA_annov23
  }
  if(!kd){
    lanno$gene_id = str_remove(lanno$gene_id,"\\.\\d*")
    manno$gene_id = str_remove(manno$gene_id,"\\.\\d*")
  }
  n1 = sum(rownames(exp) %in% manno$gene_id)
  k1 = length(n1)/nrow(exp)< 0.25 & length(n1)<5000
  n2 = sum(rownames(exp) %in% lanno$gene_id)
  k2 = length(n2)/nrow(exp)< 0.25 & length(n2)<5000
  mRNA_exp = exp[rownames(exp) %in% manno$gene_id,]
  tmp = data.frame(gene_id = rownames(exp))
  x = dplyr::inner_join(tmp,manno,by = "gene_id")
  mRNA_exp = mRNA_exp[!duplicated(x$gene_name),]
  x = x[!duplicated(x$gene_name),]
  rownames(mRNA_exp) = x$gene_name
  lnc_exp = exp[rownames(exp) %in% lanno$gene_id,]
  tmp = data.frame(gene_id = rownames(exp))
  x = dplyr::inner_join(tmp,lanno,by = "gene_id")
  lnc_exp = lnc_exp[!duplicated(x$gene_name),]
  x = x[!duplicated(x$gene_name),]
  rownames(lnc_exp) = x$gene_name
  message(paste0(nrow(mRNA_exp),
                 " of genes successfully mapping to mRNA,",
                 nrow(lnc_exp),
                 " of genes successfully mapping to lncRNA"))
  if(mrna_only){
    return(mRNA_exp)
  }else if(lncrna_only){
      return(lnc_exp)
  }else{
      return(rbind(mRNA_exp,lnc_exp))
    }
}

##' trans_array
##'
##' transform rownames for microarray expression matrix
##'
##' @param exp tcga or tcga_gtex expression set from gdc or xena
##' @param ids data.frame  with original rownames and new rownames
##' @param from colname for original rownames
##' @param to colname for new rownames
##' @return a transformed expression set with new rownames
##' @author Xiaojie Sun
##' @export
##' @examples
##' exp = matrix(1:50,nrow = 10)
##' rownames(exp) = paste0("g",1:10)
##' ids = data.frame(probe_id = paste0("g",1:10),
##'                 symbol = paste0("G",c(1:9,9)))

trans_array = function(exp,ids,from = "probe_id",
                       to = "symbol"){
  a = intersect(rownames(exp),ids[,from])
  message(paste0(length(a) ," of ",nrow(exp)," rownames matched"))
  ids = ids[!duplicated(ids[,to]),]
  exp = exp[rownames(exp) %in% ids[,from],]
  ids = ids[ids[,from]%in% rownames(exp),]
  exp = exp[ids[,from],]
  rownames(exp)=ids[,to]
  message(paste0(nrow(exp)," rownames transformed after duplicate rows removed"))
  return(exp)
}


##' sam_filter
##'
##' drop duplicated samples from the same patients
##'
##' @param exp tcga or tcga_gtex expression set from gdc or xena
##' @return a transformed expression set without duplicated samples
##' @author Xiaojie Sun
##' @export
##' @examples
##' cod2 = sam_filter(cod)

##' @seealso
##' \code{\link{simpd}};\code{\link{draw_volcano}};\code{\link{draw_venn}}

sam_filter = function(exp){
  exp = exp[,order(colnames(exp))]
  n1 = ncol(exp)
  group = make_tcga_group(exp)
  exptumor = exp[,group == "tumor"]

  expnormol = exp[,group == "normal"]
  exptumor = exptumor[,!duplicated(str_sub(colnames(exptumor),1,12))]
  expnormol = expnormol[,!duplicated(str_sub(colnames(expnormol),1,12))]

  exp = cbind(exptumor,expnormol)
  message(paste("filtered",n1-ncol(exp),"samples."))
  return(exp)
}
