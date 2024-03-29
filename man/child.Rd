% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/offspring.R
\name{child}
\alias{child}
\alias{child.tbl_tree}
\title{child}
\usage{
child(.data, .node, ...)

\method{child}{tbl_tree}(.data, .node, ...)
}
\arguments{
\item{.data}{phylo or tbl_tree object}

\item{.node}{node number}

\item{...}{additional parameters}
}
\value{
child data
}
\description{
access child data
}
\examples{
library(ape)
tree <- rtree(4)
x <- as_tibble(tree)
child(x, 4)
}
\author{
Guangchuang Yu
}
