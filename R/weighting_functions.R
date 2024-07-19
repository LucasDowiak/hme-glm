# Set of functions to obtain prior and posterior weights


#' Calculate the (unnormalized) multinomial split probabilities of an
#' arbitrary gating node given the data and the gating parameters
#' 
#' @param node The address of the gating node
#' 
#' @param gate_par_list A list (of lists) containing the parameters of the gating nodes
#' 
#' @param X The design matrix of the gating network
#' 
#' @return A matrix of un-normalized multinomial probabilities for node \code{d}
#' 
par_to_exp_gate_paths <- function(node, gate_par_list, X)
{
  pars <- gate_par_list[[node]]
  if (!is.list(pars)) {
    pars <- list(pars)
  }
  f_ <- function(beta)
  {
    return(exp(X %*% beta))
  }
  numerator <- lapply(pars, f_)
  G <- matrix(unlist(numerator), ncol=length(pars))
  G <- cbind(G, 1)
  return(G)
}


#' Calculate the multinomial split probabilities of an arbitrary gating node
#' given the data and the gating parameters
#' 
#' @param node The address of the gating node
#' 
#' @param gate_par_list A list (of lists) containing the parameters of the gating nodes
#' 
#' @param X The design matrix of the gating network
#' 
#' @return A matrix of multinomial probabilities for node \code{d}
#' 
par_to_gate_paths <- function(node, gate_par_list, X)
{
  G <- par_to_exp_gate_paths(node, gate_par_list, X)
  return(sweep(G, 1, rowSums(G), `/`))
}


#' Calculate the implied density of the observations for an arbitrary expert
#' given the data and set of regression parameters
#'
#' @param node The address of the expert node as a character
#' 
#' @param expert_type A character string specifying the type of expert
#' 
#' @param expert_par_list A list containing the parameters of the expert nodes
#' 
#' @param Y The dependent variable
#' 
#' @param X The design matrix of the expert regressions
#' 
#' @param ... Additional options to pass to the family of d*() functions
#' 
#' @return A matrix of multinomial probabilities for node \code{d}
#' 
#' @importFrom stats dnorm
#' 
par_to_expert_dens <- function(node, expert_type, expert_par_list, Y, X, ...)
{
  if (expert_type == "gaussian") {
    parm <- expert_par_list[[node]]
    beta <- parm[1:(length(parm) - 1)]
    variance <- exp(parm[length(parm)])
    mu <- X %*% beta
    return(dnorm(Y, mean=mu, sd=sqrt(variance), ...))
  } else if (expert_type == "bernoulli") {
    # if Y is a bernoulli variable {1, 0}
    g <- par_to_gate_paths(node, expert_par_list, X)
    # pmax(g[,1] * Y, g[,2] * (1 - Y))
    g[,1] * Y + g[,2] * (1 - Y)
  }
}


#' Given a list of branches retrieve the gating split probabilities 
#' 
#' Say gating node 1|0 splits into three directions: 1|0,1 - 1|0,2 - 1|0,3. This
#' function will retrieve the split probabilities of direction i from node a
#' 
#' @param branches The address of the expert node as a character
#' 
#' @param gate_prob The list containing all the current split probabilities for
#' everying gating node in the tree
#' 
#' @return A list of matrices with split probabilities
#' 
gate_path_values <- function(branches, gate_prob)
{
  if (!is(branches, "list"))
    branches <- as.list(branches)
  
  f_ <- function(b)
  {
    if (b == "0")
      return(NULL)
    branch <- unlist(last_split(b))
    node <- unlist(all_but_last_split(b))
    gate_prob[[node]][, as.integer(branch), drop=FALSE]
  }
  lapply(branches, f_)
}


#' Find the cascading split probabilities from an arbitrary gating node to one
#' of its progeny
#' 
#' @param geezer The address of the gating node closer to the root node
#' 
#' @param youngin The address of the gating node closer to the expert node
#' 
#' @param gate_prob The list containing all the current split probabilities for
#' every gating node in the tree
#' 
inter_node_paths <- function(geezer, youngin, gate_prob)
{
  # product path from node1 down to node2
  a1 <- unlist(ancestors(geezer))
  a2 <- unlist(ancestors(youngin))
  if (!geezer %in% a2) {
    stext <- "Node1 [%s] must be an ancestor of node2 [%s]."
    stop(sprintf(stext, a1, a2))
  }
  return(gate_path_values(setdiff(a2, a1), gate_prob))
}



#' Find the cumulative product path from an arbitrary gating node to one of its
#' progeny
#' 
#' @param geezer The address of the gating node closer to the root node
#' 
#' @param youngin The address of the gating node closer to the expert node
#' 
#' @param gate_prob The list containing all the current split probabilities for
#' every gating node in the tree
#' 
gate_path_product <- function(geezer, youngin, gate_prob)
{
  inp <- inter_node_paths(geezer, youngin, gate_prob)
  Reduce(`*`, inp)
}



#' Calculate the prior weight for a gating node
#' 
#' @param node An arbitrary node in the tree. Can be a gating node or an expert
#' node
#' 
#' @param gate_prob The list containing all the current split probabilities for
#' every gating node in the tree
#' 
prior_weights <- function(node, gate_prob)
{
  # prior weights go from root node down
  gate_path_product("0", node, gate_prob)
}


#' Calculate the posterior weights for a gating node
#' 
#' @param node An arbitrary node in the tree. Can be a gating node or an expert
#' node
#' 
#' @param treestr A list of the entire tree structure
#' 
#' @param gate_prob The list containing all the current split probabilities for
#' every gating node in the tree
#' 
#' @param densities The list containing the density values for all experts in
#' the tree
#' 
#' @return A matrix containing the posterior weights of the subtree starting
#' from `node`. Each column represents the posterior weight for that split
#' 
posterior_weights <- function(node, treestr, gate_prob, densities)
{
  # posterior weights go from experts up
  p <- unlist(progeny(node, treestr))
  terminals <- p[unlist(is_terminal(p, treestr))]
  childs <- unlist(children(node, treestr))
  
  Gs <- napply(terminals, function(x) gate_path_product(node, x, gate_prob))
  
  # first match terminal paths with their densities
  f_ <- function(x)
  {
    idx <- which(unlist(terminals) == x)
    Reduce(`+`, Gs[idx]) * densities[[x]]
  }
  H <- napply(names(densities), f_)
  # sum up and standardize for the appropriate number of branches
  g_ <- function(x)
  {
    Reduce(`+`, H[grepl(x, names(H))])
  }
  HH <- matrix(unlist(lapply(childs, g_)), ncol=length(childs))
  sweep(HH, 1, rowSums(HH), `/`)
}


"Calculate the joint posterior weight for a gating node (weight for the glm)

 Output:
       the cumulative product of posterior weights from the root to `node`"
#' Cumulative product of posterior weights from the root to an arbitrary node
#' 
#' @param node An arbitrary node. Can be either a gating node or an expert
#' 
#' @param post_weights The list containing all the posterior split probabilities
#' for every gating node in the tree
#' 
#' @param root_prior Weight from the branch of the parent node
#' 
joint_posterior_weight <- function(node, post_weights, root_prior)
{
  # not feasible for the root node
  lrp <- length(root_prior)
  nlp <- nrow(post_weights[[1]])
  stopifnot(lrp == 1 || lrp == nlp)
  if (node == "0") {
    if (lrp == 1)
      return(rep(1, nlp))
    else
      return(root_prior)
  }
  return(root_prior * gate_path_product("0", node, post_weights))
}


#' Calculate the likelihood contribution for a set of experts
#' 
#' @param experts A character vector of expert node names
#' 
#' @param densities he list containing the density values for all experts in
#' the tree
#' 
#' @param gate_prob The list containing all the current split probabilities for
#' every gating node in the tree
#' 
expert_lik_contr <- function(experts, densities, gate_prob)
{
  # prior weights for descendants experts
  Pi <- napply(experts, prior_weights, gate_prob)
  # prior *  P^{k}
  return(napply(experts, function(x) Pi[[x]] * densities[[x]]))
}


#' Calculates the matrix of log-likelihood contributions for all experts in the
#' tree
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
  S <- simplify2array(expert_lik_contr(expert.nodes, densities, gate_prob))
  sum(S)
}


#' Supply random initial parameter values for a gating node
#' 
#' @param node A gating node name
#' 
#' @param tree A list of the entire tree structure
#' 
#' @param n The number of parameters to provide initial values
#' 
#' @importFrom stats runif
#' 
init_gate_node_pars <- function(node, tree, n)
{
  nchilds <- length(unlist(children(node, tree)))
  limit <- 0.1
  if (nchilds == 2) {
    return(runif(n, -limit, limit))
  } else {
    return(lapply(seq_len(nchilds - 1), function(x) runif(n, -limit, limit)))
  }
}

