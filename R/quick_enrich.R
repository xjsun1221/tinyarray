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
##' @importFrom org.Hs.eg.db org.Hs.eg.db
##' @importFrom ggplot2 facet_grid
##' @importFrom ggplot2 scale_x_discrete
##' @export
##' @examples
##' gse = "GSE42872"

##' @seealso
##' \code{\link{geo_download}};\code{\link{draw_volcano}};\code{\link{draw_venn}}

quick_enrich <- function(genes,kkgo_file = "kkgo_file.Rdata"){

  if(any(is.na(as.numeric(genes)))){
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
    go <- enrichGO(genes, OrgDb = "org.Hs.eg.db", ont="all")
    save(kk,go,file = kkgo_file)
  }
  load(kkgo_file)
  kk.dot = dotplot(kk)
  go.dot = dotplot(go, split="ONTOLOGY",font.size =10,showCategory = 5)+
    facet_grid(ONTOLOGY~., scale="free") +
    scale_x_discrete(labels=function(x) str_wrap(x, width=45))
  result = list(kk = kk,go = go,kk.dot = kk.dot,go.dot = go.dot)
  return(result)
}



