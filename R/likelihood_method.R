#' Likelihood contribution for each expert
#'
#' @param experts A character vector of expert node names
#' 
#' @param densities A list containing the density values for all experts in
#' the tree
#' 
#' @param weights A list containing all the current split probabilities for
#' every gating node in the tree
#'
weighted_densities <- function(experts, densities, weights)
{
  # weights for descendants experts
  Pi <- napply(experts, root_to_node_weight, weights)
  # prior *  p^{k}
  contr <- napply(experts, function(x) Pi[[x]] * densities[[x]])
  # return as matrix
  return(simplify2array(contr, higher=FALSE))
}


#' Likelihood value for the HME
#'
#' S3 method to calculate the likelihood of a fitted HME model object. The
#' value comes from the fitted parameters of the model and not the vector of
#' historical likelihood values saved from model optimization
#'
#' @param obj An S3 object of class "hme"
#' 
#' @export
#' 
logLik.hme <- function(obj)
{
  L <- weighted_densities(obj$expert.nms, obj$list_density, obj$list_priors)
  return(sum(log(rowSums(L))))
}



#' Likelihood contributions
#' 
#' Used internally. Calculates the matrix of likelihood contributions for all
#' experts in the tree.
#' 
#' @param treestr A list of the entire tree structure
#' 
#' @param gate_prob The list containing all the current prior split probabilities
#' for every gating node in the tree
#' 
#' @param densities The list containing the density values for all experts in
#' the tree
#' 
#' @return The log-likelihood of the model
#' 
log_likelihood <- function(treestr, gate_prob, densities)
{
  expert.nodes <- treestr[unlist(is_terminal(treestr, treestr))]
  L <- weighted_densities(expert.nodes, densities, gate_prob)
  return(sum(log(rowSums(L))))
}


#' HME Selection Criterion
#'
#' Function to calculate one of several statistics that can be used for model
#' comparison.
#'
#' @param obj An S3 object of class "hme"
#' 
#' @param type One of three options: 'aic', 'bic', or 'mse'
#' 
#' @export
#' 
criterion <- function(obj, type=c("aic", "bic", "mse"))
{
  type <- match.arg(type)
  K <- obj[["no.of.pars"]]
  N <- obj[["N"]]
  L <- logLik(obj)
  
  
  if (type %in% c("aic", "bic")) {
    if (type == "aic") {
      penalty <- 2
    } else if (type == "bic") {
      penalty <- log(N)
    }
    out <- (penalty * K - 2 * L) / N
  } else if (type == "mse") {
    mses <- obj$MSE[!is.na(obj$MSE)]
    out <- mses[length(mses)]
  }
  return(out)
}
