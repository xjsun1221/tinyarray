##' quick_enrich
##'
##' do diffiencial analysis according to exprission set and group information
##'
##' @param genes a gene symbol or entrizid vector
##' @param kkgo_file Rdata filename for kegg and go result
##' @return enrichment reslut and dotplots
##' @author Xiaojie Sun
##' @importFrom clusterProfiler bitr
##' @importFrom clusterProfiler enrichKEGG
##' @importFrom clusterProfiler enrichGO
##' @importFrom clusterProfiler dotplot
##' @importFrom clusterProfiler setReadable
##' @importFrom org.Hs.eg.db org.Hs.eg.db
##' @importFrom ggplot2 facet_grid
##' @importFrom ggplot2 scale_x_discrete
##' @export
##' @examples
##' head(genes)
##' g = quick_enrich(genes)
##' names(g)
##' g[[1]][1:4,1:4]
##' g[[3]]
##' g[[4]]
##' @seealso
##' \code{\link{geo_download}};\code{\link{draw_volcano}};\code{\link{draw_venn}}

quick_enrich <- function(genes,kkgo_file = "kkgo_file.Rdata"){
  if(any(is.na(suppressWarnings(as.numeric(genes))))){
    s2e <- bitr(genes, fromType = "SYMBOL",
                toType = c( "ENTREZID"),
                OrgDb = org.Hs.eg.db::org.Hs.eg.db)
    s2e <- s2e[!duplicated(s2e$SYMBOL),]
    genes = s2e$ENTREZID
  }
  kkgo_file = kkgo_file
  if(!file.exists(kkgo_file)){
    kk <- enrichKEGG(gene         = genes,
                     organism     = 'hsa',
                     pvalueCutoff = 0.05)
    kk = setReadable(kk,OrgDb = "org.Hs.eg.db",keyType = "ENTREZID")
    go <- enrichGO(genes,
                   OrgDb = "org.Hs.eg.db",
                   ont="all",
                   readable = T)
    save(kk,go,file = kkgo_file)
  }
  load(kkgo_file)
  k1 = sum(kk@result$p.adjust <0.05)
  k2 = sum(go@result$p.adjust <0.05)
  if (k1 == 0|is.null(k1)) {
    kk.dot = "no pathway enriched"
  } else{
    kk.dot = dotplot(kk)
  }

  if (k2 == 0|is.null(k2)) {
    go.dot = "no terms enriched"
  } else{
    go.dot = dotplot(go, split="ONTOLOGY",font.size =10,showCategory = 5)+
      facet_grid(ONTOLOGY~., scale="free") +
      scale_x_discrete(labels=function(x) str_wrap(x, width=45))
  }
  result = list(kk = kk,go = go,kk.dot = kk.dot,go.dot = go.dot)
  return(result)
}



