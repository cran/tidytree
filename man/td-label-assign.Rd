% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/method-accessor.R
\name{td-label-assign}
\alias{td-label-assign}
\alias{tip.label<-}
\alias{node.label<-}
\alias{node.label<-.phylo}
\alias{node.label<-.treedata}
\alias{node.label<-.tbl_tree}
\alias{tip.label<-.phylo}
\alias{tip.label<-.treedata}
\alias{tip.label<-.tbl_tree}
\title{the tip or internal node label assign of tbl_tree phylo and treedata}
\usage{
tip.label(x) <- value

node.label(x) <- value

\method{node.label}{phylo}(x) <- value

\method{node.label}{treedata}(x) <- value

\method{node.label}{tbl_tree}(x) <- value

\method{tip.label}{phylo}(x) <- value

\method{tip.label}{treedata}(x) <- value

\method{tip.label}{tbl_tree}(x) <- value
}
\arguments{
\item{x}{object, should be one of \code{tbl_tree}, \code{phylo} or \code{treedata}}

\item{value}{character, the character vector}
}
\description{
the tip or internal node label assign of tbl_tree phylo and treedata
}
