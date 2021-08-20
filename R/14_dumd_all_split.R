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
##' @seealso
##' \code{\link{split_list}}

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

##' intersect_all
##'
##' calculate intersect  set for two or more elements
##'
##' @param ... some vectors or a list with some vectors
##' @export
##' @return vector
##' @author Xiaojie Sun
##' @examples
##' x1 = letters[1:4]
##' x2 = letters[3:6]
##' x3 = letters[3:4]
##' re =intersect_all(x1,x2,x3)
##' re2 = intersect_all(list(x1,x2,x3))
##' re3 = union_all(x1,x2,x3)
##' @seealso
##' \code{\link{union_all}}
intersect_all <- function(...){
  li = list(...)
  if(length(li)==1 & is.list(li[[1]])) li = li[[1]]
  result = li[[1]]
  for(k in 1:length(li)){
    result = intersect(result,li[[k]])
  }
  return(result)
}

##' union_all
##'
##' calculate union set for two or more elements
##'
##' @param ... some vectors or a list with some vectors
##' @export
##' @return vector
##' @author Xiaojie Sun
##' @examples
##' x1 = letters[1:4]
##' x2 = letters[3:6]
##' x3 = letters[3:4]
##' re =intersect_all(x1,x2,x3)
##' re2 = intersect_all(list(x1,x2,x3))
##' re3 = union_all(x1,x2,x3)
##' @seealso
##' \code{\link{intersect_all}}
union_all <-  function(...){
  li = list(...)
  if(length(li)==1 & is.list(li[[1]])) li = li[[1]]
  result = li[[1]]
  for(k in 1:length(li)){
    result = union(result,li[[k]])
  }
  return(result)
}


##' Split the list into sublists
##'
##' Split the list into sublists
##'
##' @param l a list with more than one element
##' @return single sublists
##' @author Xiaojie Sun
##' @export
##' @examples
##' m = list(x = 1:2,
##'          y = matrix(1:10,2),
##'          z = iris)
##' split_list(m)
##' @seealso
##' \code{\link{dumd}}

split_list = function(l){
  nam = names(l)
  if(is.null(nam)) nam = paste0("a",1:length(l))
  for(i in 1:length(l)){
    assign(nam[i],l[[i]],.GlobalEnv)
  }
}
