#' Hierarchical Mixture-of-Experts
#' 
#' An R implementation of the HME model proposed by Jordan and Jacobs
#' 
#' @param tree A character vector describing the structure of the gating network
#' 
#' @param formula A formula object in R describing the expert and network regressions.
#' The general structure is "y ~ x | z" where 'y' is the dependent variable, 'x'
#' the set of explanatory variables in the experts and 'z' is the set of variables
#' in the gating network.
#' 
#' @param data The data.frame containing the variables y, x, and z
#' 
#' @param family The glm family to use in the expert regressions
#' 
#' @param holdout An optional holdout set to track the mean squared error
#' 
#' @param init_gate_pars A list (of lists) containing starting values for the gating parameters
#' 
#' @param init_expert_pars A list containing starting values for the expert regressions
#' 
#' @param maxiter The maximum number of iterations for the EM algorithm
#' 
#' @param tol Numeric indicating the convergence tolerance for the EM algorithm
#' 
#' @param trace Control how updates of the EM algorithm are printed to the console
#' 
#' @return An object of class 'hme'
#'
#' @importFrom Formula as.Formula
#' 
#' @export
#' 
hme <- function(tree, formula, data, family=gaussian(), holdout=NULL,
                init_gate_pars=NULL, init_expert_pars=NULL,
                maxiter=100, tolerance=1e-4, trace=0)
{
  cl <- match.call()
  require(Formula)
  stopifnot(!missing(data))
  nullholdout <- is.null(holdout)
  
  # TODO: Validate Tree Structure
  tree <- sort(tree)
  hme_type <- "hme"
  expert_type <- family$family
  
  mf <- match.call(expand.dots = FALSE)
  call_ <- as.list(mf)
  m <- match(c("formula", "data"), names(mf), 0)
  form <- Formula::as.Formula(formula)
  stopifnot(length(form)[1] == 1L, length(form)[2] == 2L)
  mf <- model.frame(form, data = data)
  Y <- model.response(mf, "numeric")
  X <- model.matrix(terms(form, data = mf, rhs = 1), mf)
  Z <- model.matrix(terms(form, data = mf, rhs = 2), mf)
  N <- length(Y)
  
  if (!nullholdout) {
    mfp <- model.frame(form, data = holdout)
    Yp <- model.response(mfp, "numeric")
    Xp <- model.matrix(terms(form, data = mfp, rhs = 1), mfp)
    Zp <- model.matrix(terms(form, data = mfp, rhs = 2), mfp)
  }

  expert.nodes <- tree[unlist(is_terminal(tree, tree))]
  gate.nodes <- setdiff(tree, expert.nodes)
  expert.index <- expert_index(tree)
  
  if (is(init_gate_pars, "NULL")) {
    gate.pars <- napply(gate.nodes, init_gate_node_pars, tree, ncol(Z))
  } else {
    gate.pars <- init_gate_pars
  }
  
  if (is(init_expert_pars, "NULL")) {
    if (expert_type=="gaussian") {
      expert.pars <- bootstrap_glm(length(expert.nodes),
                                   formula=terms(form, data=mf, rhs=1),
                                   data=data,
                                   family=gaussian,
                                   model=FALSE)
      # expert.pars <- napply(expert.pars, function(x) c(x, log(var(Y))))
      names(expert.pars) <- expert.nodes
    } else {
      expert.pars <- napply(expert.nodes, function(x) runif(ncol(X) , -2, 2))
    }
  } else {
    expert.pars <- init_expert_pars
  }

  NN <- length(c(unlist(gate.pars), unlist(expert.pars)))
  
  # Run E-step on initial parameter values 
  old_priors <- napply(gate.nodes, par_to_gate_paths, gate.pars, Z)
  old_densities <- napply(expert.nodes, par_to_expert_dens, expert_type,
                           expert.pars, old_priors, Y, X)
  old_posteriors <- napply(names(gate.pars), posterior_weights,
                            tree, old_priors, old_densities)
  
  
  logL <- errorR <- matrix(NA_real_, nrow=maxiter + 1)
  parM <- matrix(NA_real_,
                 ncol=length(unlist(c(gate.pars, expert.pars))),
                 nrow=maxiter + 1)
  parM[1, ] <- unlist(c(gate.pars, expert.pars))
  logL[1, ] <- newLL <- -Inf
  
  for (ii in seq_len(maxiter)) {
    oldLL <- newLL
    mstep <- em_algo(tree, expert_type, old_posteriors,
                     Y=Y, X=X, Z=Z, exp.pars=expert.pars, gat.pars=gate.pars)
    expert.pars <- mstep[["exp.pars"]]
    gate.pars <- mstep[["gat.pars"]]
    newLL <- mstep[["loglik"]] / length(Y)
    old_posteriors <- mstep$current_posteriors
    ewmaLL <- if (is.infinite(oldLL)) {
      newLL
    } else {
      0.9 * oldLL + 0.1 * newLL
    }
    logL[ii + 1,] <- newLL
    parM[ii + 1,] <- unlist(c(gate.pars, expert.pars))
    if (!nullholdout) {
      yhat <- internal_predict(Xp, Zp, expert.nodes, gate.nodes,
                               gate.pars, expert.pars, expert_type)
      errorR[ii + 1, ] <- mean((Yp - yhat)**2)
    }
    if (trace > 0) {
      cat('\r', sprintf("Step: %d - Log-Likelihood: %f - Weighted LL: %f", ii, newLL, ewmaLL))
    } else {
      cat('\n', sprintf("Step: %d - Log-Likelihood: %f - Weighted LL: %f", ii, newLL, ewmaLL))
    }
    if (abs(oldLL - newLL) < tolerance)
      break
  }
  cat("\n")
  logL <- logL[!is.na(logL), , drop=FALSE]
  parM <- parM[apply(!is.na(parM), 1, FUN=all), ]
  errorR <- errorR[!is.na(errorR), , drop=FALSE]
  
  # Create the full score vector of theta = (omega + beta)
  gate_scores <- napply(gate.nodes, logistic_score, mstep$current_priors, mstep$current_densities, Z)
  expt_scores <- napply(expert.nodes, glm_score, expert.pars, gaussian(), mstep$current_priors,
                        mstep$current_densities, Y, X)
  scores <- do.call(cbind, c(unlist(gate_scores, recursive=FALSE), expt_scores))
  
  vcv <- sandwich_vcov(gate.nodes,
                       expert.nodes,
                       scores,
                       mstep$current_densities,
                       mstep$current_priors,
                       gate.pars,
                       expert.pars,
                       Y, X, Z, N)

  structure(list(tree=tree,
                 hme.type=hme_type,
                 gate.nodes=gate.nodes,
                 gate.pars=gate.pars,
                 gate.pars.nms=colnames(Z),
                 expert.type=expert_type,
                 expert.nms=expert.nodes,
                 expert.pars=expert.pars,
                 expert.pars.nms=c(colnames(X)), # "Dispersion"),
                 dependent.var.nme=all.vars(form)[1],
                 MSE=errorR,
                 logL=logL,
                 parM=parM,
                 list_priors=mstep$current_priors,
                 list_posteriors=mstep$current_posteriors,
                 list_density=mstep$current_densities,
                 Y=Y,
                 X=X,
                 Z=Z,
                 N=N,
                 no.of.pars=NN,
                 gate.info.matrix=NA, # deprecated
                 expert.info.matrix=NA, # deprecated
                 vcv=vcv,
                 scores=scores,
                 call_=call_),
            class="hme")
}
