##' point_cut
##'
##' calculate cut point for multiple genes
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
##' \code{\link{surv_KM}};\code{\link{surv_cox}}

point_cut <- function(exprSet_hub,meta){
  if(ncol(exprSet_hub)!=nrow(meta))stop("your exprSet_hub is not corresponds to meta")
  kp = apply(exprSet_hub, 1, function(x)length(unique(unique(x))) > 4)
  message(paste0(sum(kp)," rows with less than 5 values were ignored in cut.point calculations"))
  exprSet_hub = exprSet_hub[kp,]
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
##' calculate log_rank test p values for genes
##'
##' @inheritParams point_cut
##' @inheritParams t_choose
##' @param cut.point logical , use cut_point or not, if FALSE,use median by defult
##' @param min_gn Depending on the expression of a gene, there may be a large difference in the number of samples between the two groups, and if a smaller group of samples is less than 10 percent (default) of all, the gene will be discarded
##' @importFrom survival Surv
##' @importFrom survival survdiff
##' @export
##' @return  a vector with gene names and log_rank p value
##' @author Xiaojie Sun
##' @examples
##' surv_KM(exprSet_hub1,meta1)
##' surv_KM(exprSet_hub1,meta1,pvalue_cutoff = 1)
##' surv_KM(exprSet_hub1,meta1,cut.point = TRUE)
##' @seealso
##' \code{\link{point_cut}};\code{\link{surv_cox}}


surv_KM <- function(exprSet_hub,meta,cut.point = FALSE,pvalue_cutoff = 0.05,min_gn = 0.1){
  if(cut.point)cut_point = point_cut(exprSet_hub,meta)
  log_rank_p = c()
  for(i in 1:nrow(exprSet_hub)){
    gene = as.numeric(exprSet_hub[i,])
    if(cut.point){
      meta$group=ifelse(gene>cut_point[i],'high','low')
    }else{
      meta$group=ifelse(gene>stats::median(gene),'high','low')
    }
    data.survdiff=survival::survdiff(survival::Surv(time, event)~group,data=meta)
    log_rank_p[[i]] = 1 - stats::pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
    nn = min(table(meta$group))<= min_gn * nrow(meta)
    if(nn) log_rank_p[[i]] = 1
  }
  if(is.list(log_rank_p)) log_rank_p = unlist(log_rank_p)
  names(log_rank_p) = rownames(exprSet_hub)
  log_rank_p=sort(log_rank_p)
  tp = log_rank_p[log_rank_p<pvalue_cutoff]
  return(tp)
  message(paste0(length(tp)," of ",nrow(exprSet_hub)," genes selected"))
}

##' surv_cox
##'
##' calculate cox p values and HR for genes
##'
##' @inheritParams surv_KM
##' @param HRkeep one of "all","protect"or"risk"
##' @param continuous logical, gene expression or gene expression group
##' @importFrom survival Surv
##' @importFrom survival coxph
##' @export
##' @return a matrix with gene names ,cox p value and HR
##' @author Xiaojie Sun
##' @examples
##' surv_cox(exprSet_hub1,meta1)
##' surv_cox(exprSet_hub1,meta1,cut.point = TRUE,continuous = T)
##' surv_cox(exprSet_hub1,meta1,cut.point = TRUE,continuous = T,pvalue_cutoff = 1)
##' @seealso
##' \code{\link{point_cut}};\code{\link{surv_KM}}


surv_cox <-function(exprSet_hub,meta,cut.point = FALSE,
                    pvalue_cutoff = 0.05,HRkeep = "all",
                    continuous = FALSE,min_gn = 0.1){
  if(cut.point) cut_point = point_cut(exprSet_hub,meta)
  cox_results <-list()
  for(i in 1:nrow(exprSet_hub)){
    if(continuous) {
      gene= as.numeric(exprSet_hub[i,])
      m=survival::coxph(survival::Surv(time, event)~gene, data =  meta)
    }else{
      gene= as.numeric(exprSet_hub[i,])
      if(cut.point){
        meta$group=ifelse(gene>cut_point[i],'high','low')
      }else{
        meta$group=ifelse(gene>stats::median(gene),'high','low')
      }
      meta$group = factor(meta$group,levels = c("low","high"))
      m=survival::coxph(survival::Surv(time, event)~group, data =  meta)
      }
    beta <- stats::coef(m)
    se <- sqrt(diag(stats::vcov(m)))
    HR <- exp(beta)
    HRse <- HR * se

    #summary(m)
    tmp <- cbind(coef = beta, se = se, z = beta/se, p = 1 - stats::pchisq((beta/se)^2, 1),
                 HR = HR, HRse = HRse,
                 HRz = (HR - 1) / HRse, HRp = 1 - stats::pchisq(((HR - 1)/HRse)^2, 1),
                 HRCILL = exp(beta - stats::qnorm(.975, 0, 1) * se),
                 HRCIUL = exp(beta + stats::qnorm(.975, 0, 1) * se))
    if(continuous){
      cox_results[[i]] = tmp['gene',]
    }else{
      cox_results[[i]] = tmp['grouphigh',]
      nn = min(table(meta$group))<= min_gn * nrow(meta)
      if(nn) cox_results[[i]] = rep(NA,times = length(cox_results[[i]]))
    }
  }
  cox_results = do.call(cbind,cox_results)
  cox_results=t(cox_results)
  rownames(cox_results)= rownames(exprSet_hub)
  nn2 = apply(cox_results,1,function(x){all(is.na(x))})
  cox_results = cox_results[!nn2,]
  k = cox_results[,4]<pvalue_cutoff

  if(sum(stats::na.omit(k))==1){
    tmp = matrix(cox_results[k,],nrow = 1)
    colnames(tmp) = colnames(cox_results)
    rownames(tmp) = rownames(cox_results)[k]
    cox_results = tmp
  }else{
    cox_results=cox_results[k,]
  }

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
