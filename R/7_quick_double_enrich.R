##' quick_enrich
##'
##' do diffiencial analysis according to exprission set and group information,for human only
##'
##' @param genes a gene symbol or entrizid vector
##' @param kkgo_file Rdata filename for kegg and go result
##' @param destdir destdir to save kkgofile
##' @param species species for enrichment analysis, e.g. "human","mouse","rat"
##' @inheritParams trans_exp_new
##' @return enrichment results and dotplots
##' @author Xiaojie Sun
##' @importFrom clusterProfiler bitr
##' @importFrom clusterProfiler enrichKEGG
##' @importFrom clusterProfiler enrichGO
##' @importFrom clusterProfiler dotplot
##' @importFrom clusterProfiler setReadable
##' @importFrom ggplot2 facet_grid
##' @export
##' @examples
##' \dontrun{
##' if(requireNamespace("org.Hs.eg.db",quietly = TRUE)){
##'   head(genes)
##'   g = quick_enrich(genes,destdir = tempdir())
##'   names(g)
##'   g[[1]][1:4,1:4]
##'   g[[3]]
##'   g[[4]]
##' }else{
##'   warning("Package 'org.Hs.eg.db' needed for this function to work.
##'          Please install it by BiocManager::install('org.Hs.eg.db')",call. = FALSE)
##'   }
##' }
##' @seealso
##' \code{\link{double_enrich}}

quick_enrich <- function(genes,
                         kkgo_file = "kkgo_file.Rdata",
                         destdir = getwd(),
                         species = "human"){
  if(any(is.na(suppressWarnings(as.numeric(genes))))){
    if(species == "human"){
      if(!requireNamespace("org.Hs.eg.db",quietly = TRUE)) {
        stop("Package \"org.Hs.eg.db\" needed for this function to work.
         Please install it by BiocManager::install('org.Hs.eg.db')",call. = FALSE)
      }
      or = org.Hs.eg.db::org.Hs.eg.db
    }
    if(species == "mouse"){
      if(!requireNamespace("org.Mm.eg.db",quietly = TRUE)) {
        stop("Package \"org.Mm.eg.db\" needed for this function to work.
         Please install it by BiocManager::install('org.Mm.eg.db')",call. = FALSE)
      }
      or = org.Mm.eg.db::org.Mm.eg.db
    }
    if(species == "rat"){
      if(!requireNamespace("org.Rn.eg.db",quietly = TRUE)) {
        stop("Package \"org.Rn.eg.db\" needed for this function to work.
         Please install it by BiocManager::install('org.Rn.eg.db')",call. = FALSE)
      }
      or = org.Rn.eg.db::org.Rn.eg.db
    }
    s2e <- bitr(genes, fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = or)
    s2e <- s2e[!duplicated(s2e$SYMBOL),]
    genes = s2e$ENTREZID
  }
  f = paste0(destdir,"/",kkgo_file)
  if(!file.exists(f)){
    if(species == "human"){
      or = "org.Hs.eg.db"
      orgs = 'hsa'
    }else if(species == 'mouse'){
      or = "org.Mm.eg.db"
      orgs = 'mmu'
    }else if(species == 'rat'){
      or = "org.Rn.eg.db"
      orgs = 'rno'
    }
    kk <- enrichKEGG(gene         = genes,
                     organism     = orgs,
                     pvalueCutoff = 0.05)
    if(!is.null(kk)){kk = setReadable(kk,OrgDb = or,keyType = "ENTREZID")}
    go <- enrichGO(genes,
                   OrgDb = or,
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
      facet_grid(ONTOLOGY~., scales="free")
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
##' @param species species for enrichment analysis, e.g. "human","mouse","rat"
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
##' \dontrun{
##' if(requireNamespace("org.Hs.eg.db",quietly = TRUE)&
##'    requireNamespace("labeling",quietly = TRUE)){
##'    double_enrich(deg)
##' }else{
##'   if(!requireNamespace("org.Hs.eg.db",quietly = TRUE)) {
##'     warning("Package 'org.Hs.eg.db' needed for this function to work.
##'         Please install it by BiocManager::install('org.Hs.eg.db')",call. = FALSE)
##'   }
##'   if(!requireNamespace("labeling",quietly = TRUE)) {
##'     warning("Package 'labeling' needed for this function to work.
##'         Please install it by install.packages('labeling')",call. = FALSE)
##'   }
##' }
##' }
##' @seealso
##' \code{\link{quick_enrich}}

double_enrich <- function(deg,n = 10,color = c("#2874C5", "#f87669"), species = "human"){

  if(!requireNamespace("labeling",quietly = TRUE)) {
    stop("Package \"labeling\" needed for this function to work. Please install it byby install.packages('labeling')",call. = FALSE)
  }
  deg$change = str_to_lower(deg$change)
timestamp <- format(Sys.time(), "%Y%m%d%H%M%S")
temp_dir <- file.path(tempdir(), timestamp)
dir.create(temp_dir, recursive = TRUE)
up = quick_enrich(deg$ENTREZID[deg$change=="up"], "up.rdata", destdir = temp_dir, species = species)
down = quick_enrich(deg$ENTREZID[deg$change=="down"], "down.rdata", destdir = temp_dir, species = species)
  if(!is.null(up$kk) & !is.null(down$kk) &!is.null(up$go) &!is.null(up$go)){
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
  }else {
    warning("no pathway enriched in kegg or go,return results from quick_enrich")
    result = list(up = up,
                  down = down)
    return(result)
  }

}

utils::globalVariables(c("change","pl","Description"))


