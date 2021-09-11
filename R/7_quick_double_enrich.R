##' quick_enrich
##'
##' do diffiencial analysis according to exprission set and group information,for human only
##'
##' @param genes a gene symbol or entrizid vector
##' @param kkgo_file Rdata filename for kegg and go result
##' @param destdir destdir to save kkgofile
##' @return enrichment results and dotplots
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
##' \donttest{
##' head(genes)
##' g = quick_enrich(genes,destdir = tempdir())
##' names(g)
##' g[[1]][1:4,1:4]
##' g[[3]]
##' g[[4]]
##' }
##' @seealso
##' \code{\link{double_enrich}}

quick_enrich <- function(genes,
                         kkgo_file = "kkgo_file.Rdata",
                         destdir = getwd()){
  if(any(is.na(suppressWarnings(as.numeric(genes))))){
    s2e <- bitr(genes, fromType = "SYMBOL",
                toType = c( "ENTREZID"),
                OrgDb = org.Hs.eg.db::org.Hs.eg.db)
    s2e <- s2e[!duplicated(s2e$SYMBOL),]
    genes = s2e$ENTREZID
  }
  f = paste0(destdir,"/",kkgo_file)
  if(!file.exists(f)){
    kk <- enrichKEGG(gene         = genes,
                     organism     = 'hsa',
                     pvalueCutoff = 0.05)
    kk = setReadable(kk,OrgDb = "org.Hs.eg.db",keyType = "ENTREZID")
    go <- enrichGO(genes,
                   OrgDb = "org.Hs.eg.db",
                   ont="all",
                   readable = TRUE)
    save(kk,go,file = f)
  }
  load(f)
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
      facet_grid(ONTOLOGY~., scales="free") +
      scale_x_discrete(labels=function(x) str_wrap(x, width=45))
  }
  result = list(kk = kk,go = go,kk.dot = kk.dot,go.dot = go.dot)
  return(result)
}

##' draw enrichment bar plots for both up and down genes
##'
##' draw enrichment bar plots for both up and down genes,for human only.
##'
##' @param deg a data.frame contains at least two columns:"ENTREZID" and "change"
##' @param n how many terms will you perform for up and down genes respectively
##' @param color color for bar plot
##' @return a list with kegg and go bar plot according to up and down genes enrichment result.
##' @author Xiaojie Sun
##' @importFrom stringr str_to_lower
##' @importFrom stringr str_wrap
##' @importFrom dplyr mutate
##' @importFrom dplyr arrange
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 geom_bar
##' @importFrom ggplot2 coord_flip
##' @importFrom ggplot2 theme_light
##' @importFrom ggplot2 ylim
##' @importFrom ggplot2 scale_x_discrete
##' @importFrom ggplot2 scale_y_continuous
##' @importFrom ggplot2 theme
##' @export
##' @examples
##' \donttest{
##' double_enrich(deg)
##' }
##' @seealso
##' \code{\link{quick_enrich}}

double_enrich <- function(deg,n = 10,color = c("#2874C5", "#f87669")){

  if(!requireNamespace("labeling",quietly = TRUE)) {
    stop("Package \"labeling\" needed for this function to work. Please install it byby install.packages('labeling')",call. = FALSE)
  }
  deg$change = str_to_lower(deg$change)
  up = quick_enrich(deg$ENTREZID[deg$change=="up"],"up.rdata",destdir = tempdir())
  down = quick_enrich(deg$ENTREZID[deg$change=="down"],"down.rdata",destdir = tempdir())
  up$kk@result = mutate(up$kk@result,change = "up")
  down$kk@result = mutate(down$kk@result,change = "down")
  kk = rbind(up$kk@result[1:n,],down$kk@result[1:n,])

  up$go@result = mutate(up$go@result,change = "up")
  down$go@result = mutate(down$go@result,change = "down")
  go = rbind(up$go@result[1:n,],down$go@result[1:n,])
  ud_enrich = function(df){
    df$pl = ifelse(df$change == "up",-log10(df$p.adjust),log10(df$p.adjust))
    df = arrange(df,change,pl)
    df$Description = factor(df$Description,levels = unique(df$Description),ordered = TRUE)
    tmp = with(df, labeling::extended(range(pl)[1], range(pl)[2], m = 5))
    lm = tmp[c(1,length(tmp))]
    lm = c(floor(min(df$pl)),ceiling(max(df$pl)))
    ggplot(df, aes(x=Description, y= pl)) +
      geom_bar(stat='identity', aes(fill=change), width=.7)+
      scale_fill_manual(values = color)+
      coord_flip()+
      theme_light() +
      ylim(lm)+
      scale_x_discrete(labels=function(x) str_wrap(x, width=30))+
      scale_y_continuous(breaks = tmp,
                         labels = abs(tmp))+
      theme(
        panel.border = element_blank()
      )
  }
  result = list(kp = ud_enrich(kk),
                gp = ud_enrich(go))
  return(result)
}

utils::globalVariables(c("change","pl","Description"))


