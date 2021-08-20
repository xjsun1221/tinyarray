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
##' \code{\link{edges_to_nodes}}

interaction_to_edges = function(df,a = 1,b = 2,sep = ","){
  gs = str_split(df[,b],sep)
  edges = data.frame(a1 = rep(df[,a],times = sapply(gs,length)),
                     a2 = unlist(gs))
  edges = distinct(edges,a1,a2)
  return(edges)
}
utils::globalVariables(c("a1","a2"))

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
##' \code{\link{interaction_to_edges}}

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
