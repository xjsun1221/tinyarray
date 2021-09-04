##' plcortest
##'
##' make cor.test for given lncRNA and mRNA
##' @param lnc_exp lncRNA expression set
##' @param mRNA_exp mRNA expression set which nrow equal to lncRNA_exp
##' @param cor_cutoff cor estimate cut_off, default 0
##' @export
##' @return a list with cor.test result,names are lncRNAs, element are mRNAs
##' @author Xiaojie Sun
##' @examples
##' # to update
##' @seealso
##' \code{\link{hypertest}}

plcortest <- function(lnc_exp, mRNA_exp,cor_cutoff=0) {
  jp = list()
  for(i in 1:nrow(lnc_exp)){
    x = cbind(lnc_exp[i,],t(mRNA_exp))
    jp[[i]] = vector()
    for(j in 2: ncol(x)){
      k = stats::cor.test(x[,1],x[,j])$p.value <0.05
      k2 = stats::cor.test(x[,1],x[,j])$estimate >cor_cutoff
      if(k&k2){jp[[i]] = c(jp[[i]],colnames(x)[[j]])}
    }
  }
  names(jp) = rownames(lnc_exp)
  return(jp)
}


##' hypertest
##'
##' make hypertest for given lncRNA and mRNA common miRNAs
##' @param lnc lncRNA names
##' @param pc mRNA names
##' @param deMIR miRNA names , default NULL
##' @param lnctarget a data.frame with two column,lncRNA in the first column ,miRNA in the second column
##' @param pctarget a data.frame with two column,mRNA in the first column ,miRNA in the second column
##' @export
##' @return a data.frame with hypertest result
##' @author Xiaojie Sun
##' @examples
##' # to update
##' @seealso
##' \code{\link{plcortest}}

hypertest = function(lnc,pc,deMIR = NULL,lnctarget,pctarget){
  pcs = unique(pctarget[,1])
  pctarget = lapply(unique(pctarget[,1]), function(x){
    pctarget[,2][pctarget[,1]==x]
  })
  names(pctarget) = unique(pcs)

  lncs = unique(lnctarget[,1])
  lnctarget = lapply(unique(lnctarget[,1]), function(x){
    lnctarget[,2][lnctarget[,1]==x]
  })
  names(lnctarget) = lncs

  mir1 <- unique(unlist(lnctarget))
  mir2 <- unique(unlist(pctarget))

  mirs <- union(mir1,mir2) #合集
  popTotal <- length(mirs) #合集长度

  ceLNC <- lnc[lnc %in% names(lnctarget)] #有miRNA的lnc
  cePC <- pc[pc %in% names(pctarget)] #有miRNA的pc
  #ceMIR <- mir[mir %in% mirs]

  hyperOutput <- list()
  i <- 0
  for (lncID in ceLNC) {
    listTotal <- length(lnctarget[[lncID]]) #lnc 的miRNA数量
    for (gene in cePC) {
      i = i + 1
      ovlp <- intersect(lnctarget[[lncID]], pctarget[[gene]]) #交集

      popHits <- length(pctarget[[gene]])
      Counts <- length(ovlp)

      ovlpMIRs <- paste(ovlp, collapse = ',')
      foldEnrichment <- Counts/listTotal*popTotal/popHits
      pValue <- stats::phyper(Counts-1, popHits, popTotal-popHits,
                       listTotal, lower.tail=FALSE, log.p=FALSE)

      ceMIR <- Reduce(intersect, list(ovlp, deMIR))
      deMIRs <- paste(ceMIR, collapse = ',')
      deMIRCounts <- length(ceMIR)

      hyperOutput[[i]] <- c(lncID, gene, Counts, listTotal,
                            popHits,popTotal,foldEnrichment,pValue,ovlpMIRs,
                            deMIRCounts, deMIRs)

    }
  }

  #hyperOutput <- Reduce(rbind, hyperOutput)  ## slower
  hyperOutput <- do.call(rbind, hyperOutput)
  #hyperOutput <- rbind_list(hyperOutput) ## not test

  colnames(hyperOutput) <- c('lncRNAs','Genes','Counts','listTotal',
                             'popHits','popTotal','foldEnrichment','hyperPValue','miRNAs',
                             'deMIRCounts','deMIRs')
  hyperOutput <- as.data.frame(as.matrix(hyperOutput),
                               stringsAsFactors=FALSE)
  hyperOutput <- hyperOutput[as.numeric(hyperOutput$Counts)>0,]

  #hyperOutput$FDR <- p.adjust(as.numeric(as.character(hyperOutput$pValue)),
  #method = 'fdr')
  #hyperOutput <- hyperOutput[hyperOutput$Counts>0,]
  #hyperOutput$lncRNAs <- ensembl2symbolFun(hyperOutput$lncRNAs)
  #hyperOutput$gene <- ensembl2symbolFun(hyperOutput$gene)

  if (is.null(deMIR)) {
    hyperOutput <- hyperOutput[,! colnames(hyperOutput) %in%
                                 c('deMIRCounts','deMIRs')]
  }
  hyperOutput = hyperOutput[order(hyperOutput$hyperPValue),]
  return (hyperOutput)
  message(paste0(sum(hyperOutput$hyperPValue<0.05)," pairs p<0.05"))
}
