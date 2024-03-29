% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/offspring.R
\name{offspring}
\alias{offspring}
\alias{offspring.tbl_tree}
\title{offspring}
\usage{
offspring(.data, .node, tiponly, self_include, ...)

\method{offspring}{tbl_tree}(.data, .node, tiponly = FALSE, self_include = FALSE, ...)
}
\arguments{
\item{.data}{phylo or tbl_tree object}

\item{.node}{node number}

\item{tiponly}{whether only return tip nodes}

\item{self_include}{whether include the input node,
only applicable for tiponly = FALSE}

\item{...}{additional parameters}
}
\value{
offspring data
}
\description{
access offspring data
}
\examples{
library(ape)
tree <- rtree(4)
x <- as_tibble(tree)
offspring(x, 4)
}
\author{
Guangchuang Yu
}
