library (igraph)
library (dagitty)
library(causaleffect)
source("ProbUtils.R")

igraph_from_graphNel <- function(graphN, latNodes){
  igraph_dag <- igraph::igraph.from.graphNEL(graphN, weight = FALSE)
  for (n in latNodes) {
    adj_list <- graphN@edgeL[[n]]$edges
    if (length(adj_list) == 2) {
      igraph_dag <- igraph::add_edges(igraph_dag, c(adj_list[1], adj_list[2], adj_list[2], adj_list[1]))
      igraph_dag <- igraph::set.edge.attribute(graph = igraph_dag,
                                               name ="description",
                                               index = c(length(igraph::E(igraph_dag))-1, length(igraph::E(igraph_dag))), value = "U")
    }
  }
  for (n in latNodes){
    igraph_dag <- igraph::delete_vertices(igraph_dag, n)
  }
  return(igraph_dag)
}


############################################################
# Estimate the conditional interventional distribution     #
# of var given cond from dat and a causaleffect out object #
############################################################


getVarCondProb <- function(dat, out, vals, cdag=FALSE) {
  if (length(out$sumset) == 0  && length(out$var) > 0) {
    allVars <- out$var
    allConds <- out$cond

    if (cdag == TRUE) {
      allVars <- colnames(dat)[which(substr(colnames(dat), 1, 1) %in% out$var)]
      allConds <- colnames(dat)[which(substr(colnames(dat), 1, 1) %in% out$cond)]
    }

    #cat(paste0("returning P - ", "var: ",  allVars, "; cond: ",
    #           paste0(allConds, collapse=";"), "\n"))

    return(getCondProb(dat, allVars, allConds, vals))
  } else {
    if (length(out$sumset) > 0) { # this has to be tested before product
      #cat("computing sum...\n")
      return(computeSum(dat, out, vals, cdag))
    } else if (!is.null(out$product) && out$product) {
      #cat("computing product...\n")
      return(computeProduct(dat, out, vals, cdag))
    } else if (!is.null(out$fraction) && out$fraction) {
      return(computeFraction(dat, out, vals, cdag))
    } else {
      print("Unexpected ERROR!")
    }
  }
}


# out must have the fields num and den
computeFraction <- function(dat, out, vals, cdag=FALSE) {
  if (!is.null(out$fraction) && out$fraction) {
    num <- getVarCondProb(dat, out$num, vals, cdag)
    den <- getVarCondProb(dat, out$den, vals, cdag)

    if (den == 0)
      print("denominator has zero probability!")

    return(num/den)
  }  else {
    # ERROR
    print("Out has no product!")
  }
}


# out must have the fields, var, product, sumset
computeProduct <- function(dat, out, vals, cdag=FALSE) {
  if (!is.null(out$product) && out$product) {
    partialProd <- 1
    for (child in 1:length(out$children)) {
      #cat(paste0("computeProduct:child ", child, "\n"))
      partialProd <- partialProd * getVarCondProb(dat, out$children[[child]], vals, cdag)
    }
    return(partialProd)
  }  else {
    # ERROR
    print("Out has no product!")
  }
}

computeSum <- function(dat, out, vals, cdag=FALSE) {
  if (length(out$sumset) > 0) {

    if (cdag == FALSE) {
      sumVals <- expand.grid(rep(list(c(0,1)), length(out$sumset)))
      colnames(sumVals) <- out$sumset

      # removing fixed values of variables in the sumset
      remVars <- which(names(vals) %in% out$sumset)
    } else{
      sumsetVars <- colnames(dat)[which(substr(colnames(dat), 1, 1) %in% out$sumset)]
      sumVals <- expand.grid(rep(list(c(0,1)), length(sumsetVars)))
      colnames(sumVals) <- sumsetVars

      # removing fixed values of variables in the sumset
      remVars <- which(names(vals) %in% sumsetVars)
    }

    vals_row <- vals
    if (length(remVars) > 0) {
      vals_row <- vals[-remVars]
    }

    allSumVals <- c()
    for (i in 1:nrow(sumVals)) {
      allSumVals <- rbind(allSumVals, c(sumVals[i,, drop=F], vals_row))
    }
    sumVals <- allSumVals

    #cat("sumVals")
    #print(sumVals)

    out$sumset <- NULL
    partialSum <- 0
    for (sumRow in 1:nrow(sumVals)) {
      partialSum <- partialSum + getVarCondProb(dat, out, as.list(sumVals[sumRow,]), cdag)
    }
    return(partialSum)
  }  else {
    # ERROR
    print("Out has no sumset!")
  }
}


estimateInterDistCE <- function (dat, out, vals, cdag = FALSE) {
  return(getVarCondProb(dat, out$P, vals=vals, cdag))
}
