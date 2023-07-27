##' get_deg
##'
##' do differential analysis according to expression set and group information
##'
##' @inheritParams draw_pca
##' @inheritParams draw_volcano
##' @inheritParams trans_exp_new
##' @param entriz whether convert symbols to entriz ids
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
##' \dontrun{gse = "GSE42872"
##' geo = geo_download(gse,destdir=tempdir())
##' Group = rep(c("control","treat"),each = 3)
##' Group = factor(Group)
##' find_anno(geo$gpl)
##' ids <- AnnoProbe::idmap(geo$gpl,destdir = tempdir())
##' deg = get_deg(geo$exp,Group,ids)
##' head(deg)
##' }
##' @seealso
##' \code{\link{multi_deg}};\code{\link{get_deg_all}}

get_deg <- function(exp,
                    group_list,
                    ids,
                    logFC_cutoff=1,
                    pvalue_cutoff=0.05,
                    adjust = FALSE,
                    entriz = TRUE,
                    species = "human") {
  p3 <- is.factor(group_list)
  if(!p3) {
    group_list = factor(group_list)
    warning("group_list was covert to factor")
  }
  if(nlevels(group_list)==2){
    if(ncol(exp)!=length(group_list))stop("wrong group_list or exp")
    if(ncol(ids)!=2)stop("wrong ids pramater,it should be a data.frame with probe_id and symbol")
    colnames(ids) = c("probe_id","symbol")
    if(is.character(ids$probe_id)) ids$probe_id = as.character(ids$probe_id)

    design=stats::model.matrix(~group_list)
    fit=lmFit(exp,design)
    fit=eBayes(fit)
    deg=topTable(fit,coef=2,number = Inf)

    if("ID" %in% colnames(deg)){
      deg <- mutate(deg,probe_id=deg$ID)
    }else{
      deg <- mutate(deg,probe_id=rownames(deg))
    }
    ids = stats::na.omit(ids)
    ids = ids[!duplicated(ids$symbol),]
    deg <- inner_join(deg,ids,by="probe_id")
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

    if(entriz){
      if(species == "human"){
        if(!requireNamespace("org.Hs.eg.db",quietly = TRUE)) {
          stop("Package \"org.Hs.eg.db\" needed for this function to work.
         Please install it by BiocManger::install('org.Hs.eg.db')",call. = FALSE)
        }
        or = org.Hs.eg.db::org.Hs.eg.db
      }
      if(species == "mouse"){
        if(!requireNamespace("org.Mm.eg.db",quietly = TRUE)) {
          stop("Package \"org.Mm.eg.db\" needed for this function to work.
         Please install it by BiocManger::install('org.Mm.eg.db')",call. = FALSE)
        }
        or = org.Mm.eg.db::org.Mm.eg.db
      }
      if(species == "rat"){
        if(!requireNamespace("org.Rn.eg.db",quietly = TRUE)) {
          stop("Package \"org.Rn.eg.db\" needed for this function to work.
         Please install it by BiocManger::install('org.Rn.eg.db')",call. = FALSE)
        }
        or = org.Rn.eg.db::org.Rn.eg.db
      }

      s2e <- bitr(deg$symbol,
                  fromType = "SYMBOL",
                  toType = "ENTREZID",
                  OrgDb = or)

      deg <- inner_join(deg,s2e,by=c("symbol"="SYMBOL"))
      deg <- deg[!duplicated(deg$symbol),]
    }
  }else{
    deg <-multi_deg(exp = exp,
                    group_list = group_list,
                    ids = ids,
                    logFC_cutoff = logFC_cutoff,
                    pvalue_cutoff = pvalue_cutoff,
                    adjust = adjust,
                    entriz = entriz,
                    species = species)
  }
  return(deg)
}
