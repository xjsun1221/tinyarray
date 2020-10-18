##' draw enrichment bar plots for both up and down genes
##'
##' draw enrichment bar plots for both up and down genes
##'
##' @param deg a data.frame contains at least two columns:"ENTREZID" and "change"
##' @param n how many terms will you perform for up and down genes respectively
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
##' @importFrom labeling extended
##' @export
##' @examples
##' double_enrich(deg)
##' @seealso
##' \code{\link{draw_heatmap}};\code{\link{draw_volcano}};\code{\link{draw_venn}}

double_enrich <- function(deg,n = 10){
  deg$change = str_to_lower(deg$change)
  up = quick_enrich(deg$ENTREZID[deg$change=="up"],"up.rdata")
  down = quick_enrich(deg$ENTREZID[deg$change=="down"],"down.rdata")
  up$kk@result = mutate(up$kk@result,change = "up")
  down$kk@result = mutate(down$kk@result,change = "down")
  kk = rbind(up$kk@result[1:n,],down$kk@result[1:n,])

  up$go@result = mutate(up$go@result,change = "up")
  down$go@result = mutate(down$go@result,change = "down")
  go = rbind(up$go@result[1:n,],down$go@result[1:n,])
  ud_enrich = function(df){
    df$pl = ifelse(df$change == "up",-log10(df$p.adjust),log10(df$p.adjust))
    df = arrange(df,change,pl)
    df$Description = factor(df$Description,levels = unique(df$Description),ordered = T)
    tmp = with(df, labeling::extended(range(pl)[1], range(pl)[2], m = 5))
    lm = tmp[c(1,length(tmp))]
    lm = c(floor(min(df$pl)),ceiling(max(df$pl)))
    ggplot(df, aes(x=Description, y= pl)) +
      geom_bar(stat='identity', aes(fill=change), width=.7)+
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
