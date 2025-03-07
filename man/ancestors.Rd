% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tree_utilities.R
\name{ancestors}
\alias{ancestors}
\title{Find all ancestors of a given set of nodes}
\usage{
ancestors(d)
}
\arguments{
\item{d}{A list or character vector of nodes in the tree}
}
\value{
A list, one entry for each element of \code{d}, providing the ancestry
of each node in \code{d} back to the root node. For the root node the return
value is \code{NULL}
}
\description{
Find all ancestors of a given set of nodes
}
