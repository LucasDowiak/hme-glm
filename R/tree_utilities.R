# Set of functions to transverse the HME tree structure

#' Find children of a given set of nodes
#' 
#' @param d A list or character vector of nodes in the tree
#' 
#' @param nodes A list containing all nodes in the tree
#' 
#' @return A list, one entry for each element of \code{d}, containing the children
#' of node \code{d}. For a terminal node the return value is \code{NULL}
#' 
#' @importFrom methods is
#' 
#' @export
#' 
children <- function(d, nodes)
{
  if (!is(d, "list"))
    d <- as.list(d)
  
  f_ <- function(x)
  {
    if (is(d, "NULL"))
      return(NULL)
    o <- nodes[grep(paste0("^", x, ".[0-9]$"), nodes)]
    if (length(o) == 0) {
      return(NULL)
    } else {
      return(o)
    }
  }
  lapply(d, f_)
}


#' Find the parent of a given set of nodes
#' 
#' @param d A list or character vector of nodes in the tree
#' 
#' @param nodes A list containing all nodes in the tree
#' 
#' @return A list, one entry for each element of \code{d}, containing the parent
#' of node \code{d}. For the root node the return value is \code{NULL}
#' 
#' @export
#' 
parent <- function(d, nodes)
{
  if (!is(d, "list"))
    d <- as.list(d)
  
  f_ <- function(x)
  {
    if (x == "0") {
      return(NULL)
    } else {
      return(regmatches(x, regexpr(".*(?=\\.)", x, perl=TRUE)))
    }
  } 
  lapply(d, f_)
}


#' Find all ancestors of a given set of nodes
#' 
#' @param d A list or character vector of nodes in the tree
#' 
#' @return A list, one entry for each element of \code{d}, providing the ancestry
#' of each node in \code{d} back to the root node. For the root node the return
#' value is \code{NULL}
#' 
#' @export
#' 
ancestors <- function(d)
{
  if (!is(d, "list"))
    d <- as.list(d)
  
  fetch_history <- function(y)
  {
    f_ <- function(x, z) 
    {
      paste(c(x, z), collapse=".") 
    }
    if (y == "0") {
      return("0")
    } else {
      return(Reduce(f_, unlist(strsplit(y, "\\.")), accumulate=T))
    }
  }
  lapply(d, fetch_history)
}


#' Find the progeny of a given set of nodes
#' 
#' @param d A list or character vector of nodes in the tree
#' 
#' @param nodes A list containing all nodes in the tree
#' 
#' @return A list, one entry for each element of \code{d}, providing the ancestry
#' of each node in \code{d} back to the root node. For the root node the return
#' value is \code{NULL}
#' 
#' @export
#' 
progeny <- function(d, nodes)
{
  if (!is(d, "list"))
    d <- as.list(d)
  
  f_ <- function(z)
  {
    node_well <- c()
    look_forward <- function(y)
    {
      o <- unlist(children(y, nodes))
      if (!is(o, "NULL")) {
        node_well <<- c(o, node_well)
        look_forward(o)
      }
    }
    look_forward(z)
    if (is(node_well, "NULL")) {
      return(NULL)
    } else {
      return(sort(node_well))
    }
  }
  lapply(d, f_)
}


#' Test if a node is an expert
#' 
#' @param d A list or character vector of nodes in the tree
#' 
#' @param nodes A list containing all nodes in the tree
#' 
#' @return A list of Boolean values, one entry for each element of \code{d},
#' testing if each element of \code{d} is an expert node
#' 
#' @export
#' 
is_terminal <- function(d, nodes)
{
  if(!is(d, "list"))
    d <- as.list(d)
  
  f_ <- function(x)
  {
    bool <- FALSE
    if (is(unlist(children(x, nodes)), "NULL"))
      bool <- TRUE
    return(bool)
  }
  lapply(d, f_)
}


#' Given a node's address identify the its last split
#' 
#' @param d A list or character vector of nodes in the tree
#' 
#' @return A list
#' 
#' @export
#' 
last_split <- function(d)
{
  if(!is(d, "list"))
    d <- as.list(d)
  
  f_ <- function(x)
  {
    if (x == "0" || is(x, "NULL")) {
      return(NULL)
    } else {
      o <- unlist(strsplit(x, "\\."))
      return(o[length(o)])
    }
  }
  lapply(d, f_)
}


#' Given a node's address identify the its last split
#' 
#' @param d A list or character vector of nodes in the tree
#' 
#' @return A list
#' 
#' @export
#' 
all_but_last_split <- function(d)
{
  if(!is(d, "list"))
    d <- as.list(d)
  
  f_ <- function(x)
  {
    if (x == "0" || is(x, "NULL")) {
      return(NULL)
    } else {
      o <- unlist(strsplit(x, "\\."))
      return(paste0(o[-length(o)], collapse="."))
    }
  }
  lapply(d, f_)
}


#' Given a list of nodes return identify which ones are expert nodes
#' 
#' @param d A list or character vector of all nodes in the tree
#' 
#' @return Given a node and its address return the direction of the last split
#' 
#' @export
#' 
expert_index <- function(d)
{
  idx <- unlist(is_terminal(d, d))
  idx[idx] <- seq_len(sum(idx))
  names(idx) <- d
  return(idx)
}


#' Wrapper for lapply where x is a list of characters. Ensures the values of \code{x}
#' are passed as names to the returned list
#' 
#' @param x A list
#' 
#' @param f_ The function to apply
#' 
#' @param ... Additional variables to pass to `lapply`
#' 
#' @return A list
#' 
napply <- function(x, f_, ...)
{
  out <- lapply(x, f_, ...)
  names(out) <- x
  return(out)
}
