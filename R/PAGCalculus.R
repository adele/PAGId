#source("PAGUtils.R")

#' Applies the XUpper manipulation to the PAG, where all edges that are incoming to X are removed.
#'
#' @param amat the \code{pcalg}-compatible adjacency matrix for the PAG.
#' @param xnames the names for the variables in the set X
#'
#' @return the \code{pcalg}-compatible adjacency matrix of the
#' XUpper manipulated PAG, \eqn{P_{\overline{X}}}
#'
#' @export getXUpperManipulatedPAG
getXUpperManipulatedPAG <- function(amat, xnames) {
  amat_xupp <- amat
  for (xi_name in xnames) {
    inEdges <- which(amat_xupp[,xi_name] == 2)
    if (length(inEdges) > 0) {
      amat_xupp[inEdges, xi_name]  <- 0
      amat_xupp[xi_name, inEdges]  <- 0
    }
  }
  #plotAG(amat_xupp)
  return(amat_xupp)
}


#' Applies the XLower manipulation to the PAG, where all edges that are outgoing from X are removed.
#'
#' @param amat the \code{pcalg}-compatible adjacency matrix for the PAG.
#' @param xnames the names for the variables in the set X
#'
#' @return the \code{pcalg}-compatible adjacency matrix of the
#' XLower manipulated PAG, \eqn{P_{\underline{x}}}
#'
#' @export getXLowerManipulatedPAG
getXLowerManipulatedPAG <- function(amat, xnames) {
  amat_xlow <- amat
  for (xi_name in xnames) {
    visNodes <- getVisibleNodesFromX(amat, xi_name)
    if (length(visNodes) > 0) {
      amat_xlow[visNodes, xi_name]  <- 0
      amat_xlow[xi_name, visNodes]  <- 0
    }
  }
  #plotAG(amat_xlow)
  return(amat_xlow)
}

#' Verifies whether do-calculus' Rule 1 applies., i.e., if X and Y are
#'  m-separated by \eqn{W \cup Z} in the PAG \eqn{P_{\overline{W}}}
#'
#' @param amat the \code{pcalg}-compatible adjacency matrix for the PAG.
#' @param x the names for the variables in the set X
#' @param y the names for the variables in the set Y
#' @param z the names for the variables in the set Z
#' @param w the names for the variables in the set W
#'
#' @return a boolean indicating whether Rule 1 applies
#'
#' @export Rule1
Rule1 <- function(amat, x, y, z, w) {
  v <- colnames(amat)
  amat_upp <- getXUpperManipulatedPAG(amat, w)
  # check if all definite status paths between x and y
  # are blocked by W U Z in Pag_{\overline{W}}
  return(isMSeparated(amat_upp, x, y, c(w, z)))
}

#' Verifies whether do-calculus' Rule 2 applies., i.e., if \cr
#' X and Y are m-separated by \eqn{W \cup Z} in the PAG \eqn{P_{\overline{W}, \underline{X}}}
#'
#' @param amat the \code{pcalg}-compatible adjacency matrix for the PAG.
#' @param x the names for the variables in the set X
#' @param y the names for the variables in the set Y
#' @param z the names for the variables in the set Z
#' @param w the names for the variables in the set W
#'
#' @return a boolean indicating whether Rule 2 applies
#'
#' @export Rule2
Rule2 <- function(amat, x, y, z, w) {
  v <- colnames(amat)
  amat_upp <- getXUpperManipulatedPAG(amat, w)
  amat_lower <- getXLowerManipulatedPAG(amat_upp, x)
  # check if all definite status paths between x and y
  # are blocked by W U Z in Pag_{\overline{W}, \underline{X}}
  # print(paste0("Rule 2... x=",
  #              paste(x, collapse=","),
  #              "; y=", paste0(y, collapse=","),
  #              "; z=", paste0(z, collapse=","),
  #              "; w=", paste0(w, collapse=",")))
  return(isMSeparated(amat_lower, x, y, c(w, z)))
}

# Check if X and Y are m-separated by W U Z in Pag_{\overline{W}, \overline{X(Z)}},
#

#' Verifies whether do-calculus' Rule 3 applies., i.e., if X and Y are m-separated by
#' \eqn{W \cup Z in the PAG P_{\overline{W}, \overline{X(Z)}}, where X(Z) := X \PossAn(Z) in P_{v\setminus w}}
#'
#' @param amat the \code{pcalg}-compatible adjacency matrix for the PAG.
#' @param x the names for the variables in the set X
#' @param y the names for the variables in the set Y
#' @param z the names for the variables in the set Z
#' @param w the names for the variables in the set W
#'
#' @return a boolean indicating whether Rule 3 applies
#'
#' @export Rule3
Rule3 <- function(amat, x, y, z, w) {
  v <- colnames(amat)
  vminusw <- setdiff(v, w)
  amat_vminusw <- inducedPag(amat, vminusw)$amat
  possAnZ <- getPossAncestors(amat_vminusw, z)
  #xofZ <- setminust(x, possAnZ) # X(Z) := X \PossAn(Z)
  xofZ <- setdiff(x, possAnZ) # X(Z) := X \PossAn(Z)
  amat_upp <- getXUpperManipulatedPAG(amat, c(w, xofZ))

  # checks if all definite status paths between x and y
  # are blocked by W U Z in Pag_{\overline{W}, \overline{X(Z)}}
  return(isMSeparated(amat_upp, x, y, c(w, z)))
}

