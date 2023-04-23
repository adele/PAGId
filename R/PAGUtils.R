#' Converts a \code{pcalg} PAG into a \code{dagitty} PAG
#'
#' @param amat \code{pcalg}-compatible adjacency matrix
#'
#' @return a \code{dagitty} object for the PAG
#'
#' @import dagitty
#' @importFrom pcalg pcalg2dagitty
getDagittyPAG <- function(amat) {
  cnames <- colnames(amat)

  if (length(cnames) == 1) {
    ret_pag <- dagitty::dagitty(paste0("pag{", cnames, "}"))
  } else {
    ret_pag <- pcalg::pcalg2dagitty(amat,colnames(amat),type="pag")
  }

  return(ret_pag)
}

#' Estimates the PAG using the FCI algorithm from the \code{pcalg} R package
#' It supports datasets with all variables binary or Gaussian.
#'
#' @param dat a data-frame with the observed data.
#' @param Vnames names os the variables in the dataset, usually the column names of dat.
#' @param isbinary Boolean if TRUE, then it is expected that all variables are binary. If FALSE,
#' the it is expected that are variables follow a Gaussian distribution.
#' @param alpha the required significance level for the conditional independence tests.
#'
#' @return a \code{pcalg} object resulting from the FCI algorithm.
#'
#' @importFrom pcalg fci gaussCItest binCItest
#' @importFrom stats cor
estimatePAG <- function(dat, Vnames, isbinary, alpha) {
  ## Estimate PAG
  if (isbinary) {
    est.pag <- pcalg::fci(suffStat = list(dm = dat, adaptDF = FALSE),
                          indepTest = pcalg::binCItest, alpha = alpha, labels = Vnames, verbose = TRUE)
  } else {
    est.pag <- pcalg::fci(suffStat = list(C = stats::cor(dat), n = nrow(dat)),
                          indepTest = pcalg::gaussCItest, alpha = alpha, labels = Vnames, verbose = TRUE)
  }
  return(est.pag)
}

#' Gets the PAG induced over the variables in a set C
#'
#' @param amat the \code{pcalg}-compatible adjacency matrix of the PAG
#' @param cnames names of the variables in set C
#' @param retpag Boolean indicating whether a dagitty PAG should also be computed.
#'
#' @return a \code{pcalg}-compatible adjacency matrix for the induced PAG
inducedPag <- function (amat, cnames, retpag = TRUE) {
  ind_amat <- amat[cnames, cnames, drop = F]

  ret_pag <- NULL
  if (retpag) {
    ret_pag <- getDagittyPAG(ind_amat)
  }

  return(list("amat"=ind_amat, "pag"=ret_pag))
}

#' Gets all possible ancestors of Y in the PAG represented by the adjacency matrix amat
#'
#' @param amat the \code{pcalg}-compatible adjacency matrix of the PAG
#' @param ynames names of the variables in set Y
#'
#' @return names of all possible ancestors of Y
#'
#' @importFrom pcalg possAn
getPossAncestors <- function(amat, ynames) {
  yids <- which(colnames(amat) %in% ynames)
  an <- c()
  for (yid in yids) {
    an <- unique(c(an, pcalg::possAn(amat, x=yid, type="pag", ds=FALSE)))
  }
  an <- colnames(amat)[an]
  return(an)
}

#' Gets all possible descendants of Y in the PAG represented by the adjacency matrix amat
#'
#' @param amat the \code{pcalg}-compatible adjacency matrix of the PAG
#' @param ynames names of the variables in set Y
#'
#' @return names of all possible descendants of Y
#'
#' @importFrom pcalg possDe
getPossDescendants <- function(amat, ynames) {
  yids <- which(colnames(amat) %in% ynames)
  an <- c()
  for (yid in yids) {
    an <- unique(c(an, pcalg::possDe(amat, x=yid, type="pag", ds=FALSE)))
  }
  an <- colnames(amat)[an]
  return(an)
}



#' Gets all nodes that are connected to the set 'nodes' by
#' circle-circle edges in the PAG represented by 'amat'
#'
#' @param amat the \code{pcalg}-compatible adjacency matrix of the PAG
#' @param nodes names of a set of nodes
#' @param excludenodes names of the nodes to be excluded.
#'
#' @return nodes that are connected to the set 'nodes' by
#' circle-circle edges
getCCNodes <- function(amat, nodes, excludenodes) {
  vnames <- colnames(amat)
  ccnodes <- c()
  for (node in nodes) {
    oonodes <- names(which(amat[,node] == 1 & amat[node,] == 1))
    ccnodes <- c(ccnodes, setdiff(oonodes, excludenodes))
  }
  return(ccnodes)
}


#' Gets the nodes that are connected to X by a visible arrow out of X
#' i.e., returns \{v : x -> v is visible \}
#'
#' @param amat the \code{pcalg}-compatible adjacency matrix of the PAG
#' @param xname name of the variable X
#'
#' @return names of the nodes that are connected to X by a visible arrow out of X
#'
#' @importFrom pcalg visibleEdge
getVisibleNodesFromX <- function(amat, xname) {
  ix <- which(colnames(amat) == xname)
  potVisNodes <- which(amat[,xname] == 3) # edges with a tail in x
  visNodes <- c()
  for (node in potVisNodes) {
    if (pcalg::visibleEdge(amat, ix, node)) {
      visNodes <- c(visNodes, node)
    }
  }
  if (length(visNodes) > 0) {
    visNodes <-rownames(amat)[visNodes]
  }
  return(visNodes)
}

#' Determines if the edge v1 -> V2 is visible in the PAG
#'
#' @param amat the \code{pcalg}-compatible adjacency matrix of the PAG
#' @param v1name name of the variable V1
#' @param v2name name of the variable V2
#'
#' @return Boolean indicating if v1 -> V2 is a visible edge
#'
#' @importFrom pcalg visibleEdge
visibleEdgeByNames <- function(amat, v1name, v2name) {
  iv1 <- which(colnames(amat) == v1name)
  iv2 <- which(colnames(amat) == v2name)
  return(pcalg::visibleEdge(amat, iv1, iv2))
}

#' Gets the "invisible" possible parents of A, i.e.,
#' the set \{ V :  V o-> A or V -> A with an invisible edge\}
#'
#' @param amat_v TODO
#' @param amat the \code{pcalg}-compatible adjacency matrix of the PAG
#' @param nodeA name of the variable A
#'
#' @return names of all "invisible" possible parents of A
getInvPossParentsA <- function(amat_v, amat, nodeA) {
  invPossPar <- c()
  # nodes V such that V o-> A
  cnodes <- names(which(amat[,nodeA] == 2 & amat[nodeA,] == 1))
  invPossPar <- c(invPossPar, cnodes)

  # nodes V such that V -> A
  tnodes <- names(which(amat[,nodeA] == 2 & amat[nodeA,] == 3))
  for (tnode in tnodes) {
    visedge <- visibleEdgeByNames(amat_v, tnode, nodeA)
    if (!visedge) {
      invPossPar <- c(invPossPar, tnode)
    }
  }
  return(invPossPar)
}


#' Gets the "invisible" children of A, i.e.,
#' the set \{ V :  A o-> V or A -> V with an invisible edge\}
#'
#' @param amat_v TODO
#' @param amat the \code{pcalg}-compatible adjacency matrix of the PAG
#' @param nodeA name of the variable A
#'
#' @return names of all "invisible" possible children of A
getInvPossChildrenA <- function(amat_v, amat, nodeA) {
  invPossChild <- c()
  # nodes V such that A o-> V
  cnodes <- names(which(amat[,nodeA] == 1 & amat[nodeA,] == 2))
  invPossChild <- c(invPossChild, cnodes)

  # nodes V such that A -> V
  tnodes <- names(which(amat[,nodeA] == 3 & amat[nodeA,] == 2))
  for (tnode in tnodes) {
    if (!visibleEdgeByNames(amat_v, nodeA, tnode)) {
      invPossChild <- c(invPossChild, tnode)
    }
  }
  return(invPossChild)
}

#' Gets the possible C-Component of A, i.e.,
#'
#' @param amat_v TODO
#' @param amat the \code{pcalg}-compatible adjacency matrix of the PAG
#' @param nodesA name of the variables in the set A
#'
#' @return names of all "invisible" possible children of A
getPCComponentA <- function(amat_v, amat, nodesA) {
  pccNodesA <- c()
  for (nodeA in nodesA) {
    # nodes V such that A o-> V or A -> V with an invisible arrow
    invPossCh  <- getInvPossChildrenA(amat_v, amat, nodeA)

    newnodes <- c(nodeA, invPossCh)
    pccA <- c(newnodes)
    hhnodes <- c()
    while (length(newnodes) > 0) {
      for (node in newnodes) {
        # nodes V such that node <-> V
        hhnodes <- c(hhnodes, names(which(amat[,node] == 2 & amat[node,] == 2)))
        hhnodes <- setdiff(hhnodes, pccA)
      }
      newnodes <- hhnodes
      pccA <- c(pccA, newnodes)
    }

    possPar <- c()
    for (pccAnode in pccA) {
      # nodes V such that V o-> pccAnode or V -> pccAnode with an invisible arrow
      invPossPar <- getInvPossParentsA(amat_v, amat, pccAnode)
      possPar <- c(possPar, invPossPar)
    }
    # nodes V such that A o-o V
    ccnodes <- names(which(amat[,nodeA] == 1 & amat[nodeA,] == 1))
    pccA <- c(pccA, setdiff(c(ccnodes, possPar), pccA))
    pccNodesA <- c(pccNodesA, pccA)
  }
  return(unique(pccNodesA))
}

#' Gets the bucket of the node 'node'
#'
#' @param amat the \code{pcalg}-compatible adjacency matrix of the PAG
#' @param node name of the node
#'
#' @return the bucket of the node 'node' in the PAG represented by its adjacency matrix 'amat'
getBucket <- function(amat, node) {
  vnames <- colnames(amat)
  curnodes <- node
  bucket <- curnodes
  excludenodes <- c()
  while (length(ccnodes <- getCCNodes(amat, curnodes, excludenodes)) > 0) {
    excludenodes <-  c(excludenodes, bucket)
    bucket <- c(bucket, ccnodes)
    curnodes <- ccnodes
  }
  return(bucket)
}

#' Gets the list of all buckets in the PAG
#'
#' @param amat the \code{pcalg}-compatible adjacency matrix of the PAG
#'
#' @return a list of all buckets in the PAG represented by its adjacency matrix 'amat'
getBucketList <- function(amat) {
  vnames <- colnames(amat)
  buckets <- list()
  i <- 1
  while(length(remnodes <- setdiff(vnames, unlist(buckets))) > 0) {
    buckets[[i]] <- getBucket(amat, remnodes[1])
    i <- i+1
  }
  return(buckets)
}


#' Gets the region of a set of nodes A
#'
#' @param amat_v TODO
#' @param amat the \code{pcalg}-compatible adjacency matrix of the PAG
#' @param nodesA name of the variables in the set A
#'
#' @return the region of A
getRegion <- function(amat_v, amat, nodesA) {
  regionA <- c()
  vPCCAs <- getPCComponentA(amat_v, amat, nodesA)
  for (vPCCA in vPCCAs) {
    regionA <- unique(c(regionA, getBucket(amat, vPCCA)))
  }
  return(regionA)
}

# x is a node
getAdjNodes <- function(amat, x) {
  adjNodes <- names(which(amat[,x] != 0))
  return(adjNodes)
}

#' Determines whether X is m-separated of Y given S in the PAG represented by its adjacency matrix
#'
#' @param amat the \code{pcalg}-compatible adjacency matrix of the PAG
#' @param xnames the names for the variables in the set X
#' @param ynames the names for the variables in the set Y
#' @param snames the names for the variables in the set S
#' @param verbose a Boolean indicating whether log messages should be printed out.
#'
#' @return Boolean indicating whether X is m-separated of Y given S in the PAG
isMSeparated <- function(amat, xnames, ynames, snames, verbose=FALSE) {
  connpath <- getDefConPath(amat, xnames, ynames, snames=snames,verbose=verbose)
  return(is.null(connpath))
}

getDefConPath <- function(amat, xnames, ynames, snames, verbose=FALSE) {
  for (x in xnames) {
    for (y in ynames) {
      pathList <- list(x)
      i = 1
      while (length(pathList) > 0) {
        # there are still paths starting from X that are connecting and of definite status
        curpath <- pathList[[i]]
        len_curpath <- length(curpath)
        lastnode <- curpath[len_curpath]
        if (lastnode == y && len_curpath == 2) {
          return(curpath) # it is an edge X - Y
        }

        if (len_curpath > 2) {
          vi_name <- curpath[length(curpath)-2]
          vm_name <- curpath[length(curpath)-1]
          vj_name <- lastnode
          if (!isDefMConTriplet(amat, vi_name, vm_name, vj_name, snames)) {
            # removes the path from the list, as it will never be
            # a definite connecting path between X and Y given S
            pathList <- pathList[-i]
            next
          } else {
            # curpath starts with X and is definite connecting...
            if (lastnode == y) {
              if (verbose) {
                cat(paste0("connecting path between {",
                           paste0(xnames, collapse=","), "} and {",
                           paste0(ynames, collapse=","), "} given {",
                           paste0(snames, collapse=","), "} : ",
                           paste0(curpath, collapse=","), "\n"))
              }
              # curpath is a definite status path between X and Y is active given S
              return(curpath)
            }
          }
        }
        # curpath has less than three nodes or is active given S
        # adding neighbors to see if it is a path to Y
        lastnode_i <- which(colnames(amat) == lastnode)
        adjnodes <- setdiff(getAdjNodes(amat, lastnode), curpath)
        n_adjn <- length(adjnodes)
        if (n_adjn == 0) {
          pathList <- pathList[-i] # path does not reach y
        } else {
          # checking if Y is one of the adj nodes
          if (y %in% adjnodes) {
            # replaces the current path with the one that has additionally the node Y
            pathList[[i]] <- c(curpath, y)
            adjnodes <- setdiff(adjnodes, y)
            n_adjn <- n_adjn - 1
          } else {
            pathList <- pathList[-i]
          }
          if (n_adjn > 0) {
            temp <- matrix(rep(curpath, n_adjn), ncol=n_adjn)
            temp <- rbind(temp, adjnodes)
            pathList <- append(pathList, split(temp, rep(1:ncol(temp), each = nrow(temp))))
          }
        }
      }
    }
  }
  # All paths between xnames and ynames are m-separated by snames
  return(NULL)
}


# triplet is an array of three nodes
isCollider <- function(amat, triplet) {
  if (length(triplet) != 3) {
    print("A triplet must have exactly three nodes")
    return(NULL)
  }
  vi <- triplet[1]
  vm <- triplet[2]
  vj <- triplet[3]
  return(amat[vi, vm] == 2 && amat[vj, vm] == 2)
}

# triplet is an array of three nodes
isDefiniteNonCollider <- function(amat, triplet) {
  if (length(triplet) != 3) {
    print("A triplet must have exactly three nodes")
    return(NULL)
  }
  vi <- triplet[1]
  vm <- triplet[2]
  vj <- triplet[3]
  return ( (amat[vi, vm] == 3 || amat[vj, vm] == 3) ||  # vm has a tail
             (amat[vi, vm] == 1 && amat[vj, vm] == 1 && amat[vi, vj] == 0) ) # vi -o vm o- vj
}

# apath is a sequence of nodes connected by an edge
isDefiniteMConnecting <- function(amat, apath, Z) {
  n <- length(apath)
  if (n < 3) {
    return(TRUE)
  }
  vnames <- colnames(amat)
  for (i in 1:(n-2)) {
    vi_name <- apath[i]
    vm_name <- apath[i+1]
    vj_name <- apath[i+2]
    if (!isDefMConTriplet(amat, vi_name, vm_name, vj_name, Z)) {
      # triplet is blocked or has a non-definite status
      return(FALSE)
    }
  }
  return(TRUE)
}

isDefMConTriplet <- function(amat, vi_name, vm_name, vj_name, Z) {
  vnames <- colnames(amat)
  connecting <- TRUE
  triplet <- c(vi_name, vm_name, vj_name)
  if (isCollider(amat, triplet)) {
    vm_id <- which(vnames == triplet[2])
    devm <- vnames[pcalg::searchAM(amat, vm_id, type="de")]
    if (length(which(devm %in% Z)) == 0) {
      # no descendant of this collider is in Z
      connecting <- FALSE
    }
  } else if (isDefiniteNonCollider(amat, triplet)) {
    if (length(which(Z == vm_name)) > 0) {
      # this non-collider is in Z
      connecting <- FALSE
    }
  } else {
    connecting <- FALSE # triplet has not a definite status
  }
  return(connecting)
}
