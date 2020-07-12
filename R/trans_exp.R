##' trans_exp
##'
##' transform rownames of tcga or tcga_gtex expression set from gdc or xena,from ensamble id to gene symbol
##'
##' @param exp tcga or tcga_gtex expression set from gdc or xena
##' @param mrna_only only keep mrna rows in result
##' @param lncrna_only only keep lncrna rows in result
##' @return a transformed expression set with symbol
##' @author Xiaojie Sun
##' @importFrom stringr str_detect
##' @importFrom dplyr inner_join
##' @export
##' @examples
##' exp = matrix(rnorm(1000),ncol = 10)
##' rownames(exp) = sample(mRNA_annov23$gene_id,100)
##' colnames(exp) = c(paste0("TCGA",1:5),paste0("GTEX",1:5))
##' k  = trans_exp(exp)

##' @seealso
##' \code{\link{simpd}};\code{\link{draw_volcano}};\code{\link{draw_venn}}

trans_exp = function(exp,mrna_only = F,lncrna_only = F){
  k00 = any(str_detect(colnames(exp),"TCGA"))
  if(!k00)warning("this expression set probably not from TCGA,please ensure it")
  k0 = any(str_detect(colnames(exp),"GTEX"))
  if(!k0){
    lanno = lnc_anno
    manno = mRNA_anno
  }else{
    lanno = lnc_annov23
    manno = mRNA_annov23
    }
  n1 = sum(rownames(exp) %in% manno$gene_id)
  k1 = length(n1)/nrow(exp)< 0.25 & length(n1)<5000
  n2 = sum(rownames(exp) %in% lanno$gene_id)
  k2 = length(n2)/nrow(exp)< 0.25 & length(n2)<5000
  message(paste0(n1,
                 " of genes successfully mapping to mRNA,",
                 n2,
                 " of genes successfully mapping to lncRNA"))
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
  if(mrna_only){
    return(mRNA_exp)
  }else if(lncrna_only){
      return(lnc_exp)
  }else{
      return(rbind(mRNA_exp,lnc_exp))
    }
  }
