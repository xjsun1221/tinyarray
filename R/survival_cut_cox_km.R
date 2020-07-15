##' point_cut
##'
##' caculate cut point for mutiple genes
##'
##' @param exprSet_hub a tumor expression set for hubgenes
##' @param meta meta data corresponds to expression set
##' @importFrom survminer surv_cutpoint
##' @export
##' @return a vector with cutpoint for genes
##' @author Xiaojie Sun
##' @examples
##' point_cut(exprSet_hub1,meta1)
##' @seealso
##' \code{\link{geo_download}};\code{\link{draw_volcano}};\code{\link{draw_venn}}

point_cut <- function(exprSet_hub,meta){
  if(ncol(exprSet_hub)!=nrow(meta))stop("your exprSet_hub is not corresponds to meta")
  dat = cbind(t(exprSet_hub),meta)
  cut.point = c()
  for(i in 1:nrow(exprSet_hub)){
    cut = survminer::surv_cutpoint(
      dat,
      time = "time",
      event = "event",
      variables = rownames(exprSet_hub)[i]
    )
    cut.point[[i]] = cut[["cutpoint"]][1,1]
  }
  names(cut.point) = rownames(exprSet_hub)
  return(cut.point)
}

##' surv_KM
##'
##' caculate log_rank test p values for genes
##'
##' @inheritParams point_cut
##' @inheritParams t_choose
##' @param cut.point logical , use cut_point or not, if FALSE,use median by defult
##' @importFrom survival Surv
##' @importFrom survival survdiff
##' @export
##' @return a vector with gene names and log_rank p value
##' @author Xiaojie Sun
##' @examples
##' surv_KM(exprSet_hub1,meta1)
##' surv_KM(exprSet_hub1,meta1,cut.point = T)
##' @seealso
##' \code{\link{geo_download}};\code{\link{draw_volcano}};\code{\link{draw_venn}}


surv_KM <- function(exprSet_hub,meta,cut.point = F,pvalue_cutoff = 0.05){
  cut_point = point_cut(exprSet_hub,meta)
  log_rank_p = c()
  for(i in 1:nrow(exprSet_hub)){
    gene = as.numeric(exprSet_hub[i,])
    if(cut.point){
      meta$group=ifelse(gene>cut_point[i],'high','low')
    }else{
      meta$group=ifelse(gene>median(gene),'high','low')
    }
    data.survdiff=survival::survdiff(survival::Surv(time, event)~group,data=meta)
    log_rank_p[[i]] = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  }
  names(log_rank_p) = rownames(exprSet_hub)
  log_rank_p=sort(log_rank_p)
  tp = log_rank_p[log_rank_p<pvalue_cutoff]
  return(tp)
  message(paste0(length(tp)," of ",nrow(exprSet_hub)," genes selected"))
}

##' surv_cox
##'
##' caculate cox p values and HR for genes
##'
##' @inheritParams surv_cox
##' @param HRkeep one of "all","protect"or"risk"
##' @importFrom survival Surv
##' @importFrom survival coxph
##' @export
##' @return a matrix with gene names ,cox p value and HR
##' @author Xiaojie Sun
##' @examples
##' surv_cox(exprSet_hub1,meta1,cut.point = T,HRkeep = "all")
##' @seealso
##' \code{\link{geo_download}};\code{\link{draw_volcano}};\code{\link{draw_venn}}


surv_cox <-function(exprSet_hub,meta,cut.point = F,pvalue_cutoff = 0.05,HRkeep = "all"){
  cut_point = point_cut(exprSet_hub,meta)
  cox_results <-list()
  for(i in 1:nrow(exprSet_hub)){
    #i = 1
    gene= as.numeric(exprSet_hub[i,])
    if(cut.point){
      meta$group=ifelse(gene>cut_point[i],'high','low')
    }else{
      meta$group=ifelse(gene>median(gene),'high','low')
    }
    meta$group = factor(meta$group,levels = c("low","high"))
    m=survival::coxph(survival::Surv(time, event) ~　group, data =  meta)

    beta <- coef(m)
    se <- sqrt(diag(vcov(m)))
    HR <- exp(beta)
    HRse <- HR * se

    #summary(m)
    tmp <- round(cbind(coef = beta, se = se, z = beta/se, p = 1 - pchisq((beta/se)^2, 1),
                       HR = HR, HRse = HRse,
                       HRz = (HR - 1) / HRse, HRp = 1 - pchisq(((HR - 1)/HRse)^2, 1),
                       HRCILL = exp(beta - qnorm(.975, 0, 1) * se),
                       HRCIUL = exp(beta + qnorm(.975, 0, 1) * se)), 3)
    cox_results[[i]] = tmp['grouphigh',]
  }
  cox_results = do.call(cbind,cox_results)
  cox_results=t(cox_results)
  rownames(cox_results)= rownames(exprSet_hub)
  cox_results=cox_results[cox_results[,4]<pvalue_cutoff,]
  if(HRkeep == "all"){
    return(cox_results)
    message(paste0(nrow(cox_results)," of ",nrow(exprSet_hub)," genes selected"))
  }else if(HRkeep == "protect"){
    k = cox_results[,5]<1
    return(cox_results[k])
    message(paste0(sum(k)," of ",nrow(exprSet_hub)," genes selected"))
  }else if(HRkeep == "risk"){
    k = cox_results[,5]>1
    return(cox_results[k])
    message(paste0(sum(k)," of ",nrow(exprSet_hub)," genes selected"))
  }else if(!(HRkeep %in% c("all","protect","risk"))){
    stop('HRkeep should be one of "all","protect"or"risk"')
  }
}
