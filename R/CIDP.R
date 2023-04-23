#' @title Algorithm for identifying conditional causal effects from PAGs
#'
#' @description Decides whether the interventional distribution P(y|do(x),z) is
#' identifiable or not given a PAG and the observational distribution P(V).
#'
#' @param amat the pcalg adjacency matrix for the PAG.
#' @param x the names for the variables in the treatment set X
#' @param y the names for the variables in the outcome set Y
#' @param z the names for the variables in the covariate set Z
#' @param verbose a Boolean indicating whether log messages should be printed out.
#'
#' @return a named list where
#'   'id'    - a boolean indicating whether P(y|do(x)) is identifiable or not;
#'   'Qexpr' - list with every required operation in the identification formula as a string;
#'   'Qop'   - list with every required operation in the identification formula;
#'   'query' - a string representing P(y|do(x)). It can be used to access the high-level
#'    identification expression for the query (e.g.,. Qexpr[[query]] and Qop[[query]])
#'
#' @seealso \code{\link{IDP}}
#'
#' @export CIDP
CIDP <- function(amat, x, y, z, verbose=TRUE) {
  if (is.null(z)) {
    return(IDP(amat, x, y))
  }

  checkCondLine2 <- function(buckets, d) {
    for (bi in buckets) {
      if (length(intersect(bi, d)) > 0 && !all(bi %in% d)) {
         return(bi)
      }
    }
    return(NULL)
  }

  checkCondLine9 <- function(x, y, z, zpartition) {
    if (length(zpartition) > 0) {
      for (i in 1:length(zpartition)) {
        zi <- zpartition[[i]]
        if (Rule2(amat, zi, y, setdiff(z, zi), x)) {
          return(i)
        }
      }
    }
    return(NULL)
  }

  v <- colnames(amat)

  intervset <- x
  intervset_str <- "obs"
  Pinterv <- "P"
  if( length(intervset) > 0) {
    intervset_str <- paste0(intervset, collapse=",")
    Pinterv <- paste0("P_{", intervset_str, "}")
  }
  query <- paste0(Pinterv, "(", paste0(y, collapse=","), " | ", paste0(z, collapse=","), ")")

  vminusx <- setdiff(v, x)
  amat_vminusx <- inducedPag(amat, vminusx)$amat
  # Computing PossAn(Y U Z) in P_vminusx
  d <- getPossAncestors(amat_vminusx, c(y,z))
  # Buckets in P
  buckets <- getBucketList(amat)

  while(!is.null(bi <- checkCondLine2(buckets, d) )) {
    # bi \cap d != 0 and bi is not contained in d
    if (verbose)
      print(paste0("bi={", paste0(bi, collapse = ","), "} satisfies Cond Line 2"))
    xprime <- intersect(bi, x)
    w <- setdiff(x, xprime)
    if (Rule2(amat, xprime, y, z, w)) {
      x <- w
      z <- union(z, xprime)
      vminusx <- setdiff(v, x)
      amat_vminusx <- inducedPag(amat, vminusx)$amat
      d <- getPossAncestors(amat_vminusx, c(y, z))
    } else {
      # throw FAIL
      if (verbose)
        print(paste0("FAIL in Line 8 for B={",
                  paste0(bi, collapse=","), "} and D={", paste0(d, collapse=","), "}"))
      return(list(id=FALSE, query=query))
    }
  }

  zpartition <- list()
  for (bi in buckets) {
    inter <- c(intersect(bi, z))
    if (length(inter) > 0)
      zpartition <- append(zpartition, list(inter))
  }

  while(!is.null(zpid <- checkCondLine9(x, y, z, zpartition) )) {
    zi <- zpartition[[zpid]]
    if (verbose)
      print(paste0("zi={", paste0(zi, collapse = ","), "} satisfies Cond Line 9"))
    x <- union(x, zi)
    z <- setdiff(z, zi)
    zpartition <- zpartition[-zpid]
  }

  ret <- IDP(amat, x, c(y, z)) # P_x(y,z)

  if (ret$id) {
    # P_x(y|z) = ret$query / \sum_y ret$query
    den <- paste0("\\sum_{", paste0(y, collapse=","), "}", ret$query)
    ret$Qexpr[[query]] <- paste0("\\frac{", ret$query,"}{", den, "}")
    ret$Qop[[query]] <- list(type="frac_cond", param=list(den.sumset=y, prob=ret$query))
    return(list(id=TRUE, query=query, Qop=ret$Qop, Qexpr=ret$Qexpr))
  } else {
    return(list(id=FALSE, query=query))
  }
}



#source("PAGUtils.R")
#source("PAGCalculus.R")

# Returns P_x(y|z) or FAIL
# Input:
# amat: a pcalg adj matrix corresponding to apag over V
