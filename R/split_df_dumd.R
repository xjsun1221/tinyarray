##' split_df
##'
##' split one column by sep paramter,replicate anothor column to return a two-column df
##'
##' @param df data.frame with at least two column
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
##' split_df(df)
##' @seealso
##' \code{\link{geo_download}};\code{\link{draw_volcano}};\code{\link{draw_venn}}

split_df = function(df,a = 1,b = 2,sep = ","){
  gs = str_split(df[,b],sep)
  sdf = data.frame(a1 = rep(df[,a],times = sapply(gs,length)),
                   a2 = unlist(gs))
  sdf = distinct(a1,a2)
  return(sdf)
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

