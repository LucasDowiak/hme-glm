% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/weighting_functions.R
\name{posterior_weights}
\alias{posterior_weights}
\title{Calculate the posterior weights for a gating node}
\usage{
posterior_weights(node, treestr, gate_prob, densities)
}
\arguments{
\item{node}{An arbitrary node in the tree. Can be a gating node or an expert
node}

\item{treestr}{A list of the entire tree structure}

\item{gate_prob}{The list containing all the current split probabilities for
every gating node in the tree}

\item{densities}{The list containing the density values for all experts in
the tree}
}
\value{
A matrix containing the posterior weights of the subtree starting
from `node`. Each column represents the posterior weight for that split
}
\description{
Calculate the posterior weights for a gating node
}
