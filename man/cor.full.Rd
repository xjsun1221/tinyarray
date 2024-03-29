% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/9_cor_t_choose.R
\name{cor.full}
\alias{cor.full}
\title{cor.test for all variables}
\usage{
cor.full(x, drop = min(x) - 0.001, min.obs = 10)
}
\arguments{
\item{x}{A numeric matrix or data.frame}

\item{drop}{drop values}

\item{min.obs}{minimum number of observations  after dropping}
}
\value{
a data.frame with cor.test p.value and estimate
}
\description{
cor.test for all variables(each two columns)
}
\examples{
x = iris[,-5]
cor.full(x)
}
\seealso{
\code{\link{cor.one}}
}
\author{
Xiaojie Sun
}
