##' interaction_to_edges
##'
##' split interactions by sep paramter,return edges data.frame
##'
##' @param df interactions data.frame
##' @param a column to replicate
##' @param b column to split
##' @param sep a character string to separate b column
##' @importFrom stringr str_split
##' @importFrom dplyr distinct
##' @export
##' @return a new df with two column ,one interaction by one rows
##' @author Xiaojie Sun
##' @examples
##' df = data.frame(a = c("gene1","gene2","gene3"),
##' b = c("d,f,a,b",
##' "c,e,g",
##' "a,b,d"))
##' interaction_to_edges(df)
##' @seealso
##' \code{\link{geo_download}};\code{\link{draw_volcano}};\code{\link{draw_venn}}

interaction_to_edges = function(df,a = 1,b = 2,sep = ","){
  gs = str_split(df[,b],sep)
  edges = data.frame(a1 = rep(df[,a],times = sapply(gs,length)),
                   a2 = unlist(gs))
  edges = distinct(edges,a1,a2)
  return(edges)
}

##' edges_to_nodes
##'
##' get nodes from edges
##'
##' @param edges data.frame
##' @importFrom stringr str_split
##' @importFrom dplyr distinct
##' @export
##' @return nodes data.frame
##' @author Xiaojie Sun
##' @examples
##' df = data.frame(a = c("gene1","gene2","gene3"),
##' b = c("d,f,a,b",
##' "c,e,g",
##' "a,b,d"))
##' edges = interaction_to_edges(df)
##' nodes = edges_to_nodes(edges)
##' @seealso
##' \code{\link{geo_download}};\code{\link{draw_volcano}};\code{\link{draw_venn}}

edges_to_nodes = function(edges){
  if(!is.null(colnames(edges))){
    m = colnames(edges)
  }else{
    m = c("A1","A2")
  }
  a = unique(edges[,1])
  b = unique(edges[,2])
  nodes = data.frame(gene = c(a,b),
                     type = c(rep(m[1],times = length(a)),
                              rep(m[2],times = length(b))))
  return(nodes)
}


##' count unique values in every colunms for data.frame
##'
##' in geo analysis,this function can help you simplify pdata, delete columns with unique values,which can't be used as group vector
##' @param x A data.frame.
##' @return The simple data.frame of columns unique values count in \code{x}
##' @importFrom dplyr arrange
##' @importFrom dplyr desc
##' @importFrom tibble tibble
##' @importFrom dplyr %>%
##' @export
##' @examples
##' dumd(iris)
##' data(ToothGrowth)
##' x = ToothGrowth
##' dumd(ToothGrowth)
##' @section just :
##' See what are you doing

dumd <- function(x){
  colname <- vector("character")
  count <- vector("integer")
  for(i in 1:ncol(x)){
    colname[i] = colnames(x)[[i]]
    count[i]=nrow(x[!duplicated(x[,i]),])
  }
  df <- tibble(colname,count) %>%
    arrange(desc(count))
  print(df)
}

