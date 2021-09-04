##' t_choose
##'
##' choose differential expressed genes by simple t.test
##'
##' @inheritParams get_deg
##' @param genes a vector with some genes
##' @param up_only keep up genes in the result only
##' @param down_only keep down genes in the result only
##' @param pvalue_cutoff p value cut off ,0.05 by defult
##' @export
##' @return a vector with differential expressed genes
##' @author Xiaojie Sun
##' @examples
##' exp = matrix(rnorm(1000),ncol = 10)
##' rownames(exp) = sample(mRNA_annov23$gene_id,100)
##' colnames(exp) = c(paste0("TCGA",1:5),paste0("GTEX",1:5))
##' exp2  = trans_exp(exp)
##' exp2[,1:5] = exp2[,1:5]+10
##' group_list = rep(c("A","B"),each = 5)
##' genes = sample(rownames(exp2),3)
##' t_choose(genes,exp2,group_list)

t_choose <- function(genes,exp,group_list,up_only = FALSE,down_only = FALSE,pvalue_cutoff = 0.05){
  if(up_only&down_only)stop("please change neither up_only or down_only to FALSE")
  genes = genes[genes %in% rownames(exp)]
  exp_small = exp[genes,]
  dat = data.frame(t(exp_small),check.names = FALSE)
  dat$group_list = group_list
  p_v <- sapply(1:(ncol(dat)-1), function(i){
    stats::t.test(dat[,i] ~group_list)$p.value
  })
  names(p_v) = colnames(dat)[-ncol(dat)]

  exp_genes = names(p_v[p_v < pvalue_cutoff])

  if(up_only){
    es_up <- sapply(1:(ncol(dat)-1), function(i){
      tmp = stats::t.test(dat[,i] ~group_list)
      k = tmp$estimate[2]-tmp$estimate[1] >0
      return(k)
    })
    up_genes = names(p_v)[p_v < pvalue_cutoff & es_up]
    return(up_genes)
  }else if(down_only){
    es_down <- sapply(1:(ncol(dat)-1), function(i){
      tmp = stats::t.test(dat[,i] ~group_list)
      k = tmp$estimate[2]-tmp$estimate[1] <0
      return(k)
    })
    down_genes = names(p_v)[p_v <pvalue_cutoff & es_down]
    return(down_genes)
  }else{
    return(exp_genes)
  }
}


##' cor.test for all variables
##'
##' cor.test for all variables(each two columns)
##'
##' @param x A numeric matrix or data.frame
##' @return a data.frame with cor.test p.value and estimate
##' @author Xiaojie Sun
##' @export
##' @examples
##' x = iris[,-5]
##' cor.full(x)
##' @seealso
##' \code{\link{cor.one}}

cor.full <- function(x){
  ss = list()
  p = list()
  ss1 = utils::combn(colnames(x),2)
  ss2 = apply(ss1, 2, paste,collapse =":")

  for(i in (1:ncol(ss1))){
    bt = x[,ss1[1,i]]
    kt = x[,ss1[2,i]]
    cot = stats::cor.test(bt,kt)
    p[[i]] = c(cot$p.value,cot$estimate)
    names(p[[i]]) = c("p.value","cor")
  }
  re = do.call(cbind,p)
  colnames(re) = apply(ss1, 2, paste,collapse =":")
  return(as.data.frame(t(re)))
}



##' cor.test for one variable with all variables
##'
##' cor.test for all variables(each two columns)
##'
##' @param x A numeric matrix or data.frame
##' @param var your chosen variable,only one.
##' @return A data.frame with cor.test p.value and estimate
##' @author Xiaojie Sun
##' @export
##' @examples
##' x = iris[,-5]
##' cor.one(x,"Sepal.Width")
##' @seealso
##' \code{\link{cor.full}}

cor.one <- function(x,var){
  if(!(var %in% colnames(x))) stop(paste0(var," is not a colname of ",x,",please check it."))
  if(!all(!duplicated(colnames(x)))) stop("unique colnames is required")
  ss = list()
  p = list()
  ss1 = matrix(c(rep(var,times = (ncol(x)-1)),
                 setdiff(colnames(x),var)),
               nrow = 2,byrow = TRUE)
  ss2 = setdiff(colnames(x),var)

  for(i in (1:ncol(ss1))){
    bt = x[,ss1[1,i]]
    kt = x[,ss1[2,i]]
    cot = stats::cor.test(bt,kt)
    p[[i]] = c(cot$p.value,cot$estimate)
    names(p[[i]]) = c("p.value","cor")
  }
  re = do.call(cbind,p)
  colnames(re) = ss2
  return(as.data.frame(t(re)))
}

