##' t_choose
##'
##' choose diffrencial expressed genes by simple t.test
##'
##' @inheritParams get_deg
##' @param genes a vector with some genes
##' @param up_only keep up genes in the result only
##' @param down_only keep down genes in the result only
##' @export
##' @return a vector with diffrencial expressed genes
##' @author Xiaojie Sun
##' @examples
##' exp = matrix(rnorm(1000),ncol = 10)
##' rownames(exp) = sample(mRNA_annov23$gene_id,100)
##' colnames(exp) = c(paste0("TCGA",1:5),paste0("GTEX",1:5))
##' k  = trans_exp(exp)
##' k[,1:5] = k[,1:5]+10
##' exp = k
##' group_list = rep(c("A","B"),each = 5)
##' genes = sample(rownames(k),3)
##' t_choose(genes,k,group_list)
##' @seealso
##' \code{\link{geo_download}};\code{\link{draw_volcano}};\code{\link{draw_venn}}


t_choose <- function(genes,exp,group_list,up_only = F,down_only = F){
  table(genes %in% rownames(exp))
  exp_small = exp[rownames(exp) %in% genes,]
  dat = data.frame(t(exp_small))
  dat$group_list = group_list
  p_v <- sapply(1:(ncol(dat)-1), function(i){
    t.test(dat[,i] ~group_list)$p.value
  })
  names(p_v) = colnames(dat)[-ncol(dat)]
  table(p_v<0.05)
  exp_genes = names(p_v[p_v < 0.05])
  if(up_only){
    es_up <- sapply(1:(ncol(dat)-1), function(i){
      tmp = t.test(dat[,i] ~group_list)
      k = tmp$estimate[2]-tmp$estimate[1] >0
      return(k)
    })
    up_genes = genes[p_v <0.05 & es_up]
    return(up_genes)
  }else if(down_only){
    es_down <- sapply(1:(ncol(dat)-1), function(i){
      tmp = t.test(dat[,i] ~group_list)
      k = tmp$estimate[2]-tmp$estimate[1] <0
      return(k)
    })
    down_genes = genes[p_v <0.05 & es_down]
    return(down_genes)
  }else{
    return(exp_genes)
  }
}
