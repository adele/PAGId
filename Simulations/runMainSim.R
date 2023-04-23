# install.packages("BiocManager")
# BiocManager::install("RBGL")
# BiocManager::install("graph")
# BiocManager::install("Rgraphviz")

# install.packages("igraph", dependencies=TRUE)
# install.packages("dagitty", dependencies=TRUE)
# install.packages("pcalg", dependencies=TRUE)
# install.packages("bnlearn", dependencies=TRUE)
# install.packages("causaleffect", dependencies=TRUE)
# install.packages("PAGId", dependencies=TRUE)

library(pcalg)
library(dagitty)
library(igraph)
library(causaleffect)
library(bnlearn)
library(PAGId)
source("PAGInstances.R")
source("ProbUtils.R")
source("EstimateCausalEffect.R")

library(doParallel)
registerDoParallel(cl <- makeCluster(8)) # running in 8 cores

# returns all combinations of adag$x and adag$z with y=1
getDefaultValsGrid <- function(adag) {
  vnames <- c(adag$x, adag$z)
  vals_grid <- expand.grid(comprehenr::to_list(for (var in vnames) c(0,1)))
  colnames(vals_grid) <- vnames
  vals_grid[[adag$y]] <- 1
  return(vals_grid)
}


runTest <- function(type, N) {
  vals_grid <- NULL
  if (type == "g1") {
    cat("Testing Example 1, P_{w,x1,x2}(y|z) is identifiable\n\n")
    adag <- getDAG1() # corresponding PAG is the same as getPAG1()
  } else if (type == "g2n") {
    cat("Testing Example 2n, P_{a,f}(y|b,e) is identifiable\n\n")
    adag <- getDAG2n()
  } else if (type == "g2") {
    cat("Testing Example 2, P_c(y|b, e) is identifiable\n\n")
    adag <- getDAG2() # corresponding PAG is the same as getPAG2()
  } else if (type == "2bR42") {
    cat("Testing Example 2b in R42, P_{x1, x2}(y1, y2, y3, y4, y5) is identifiable\n\n")
    adag <- getDAG2bR42()
    vnames <- c(adag$y, adag$x) # testing all combinations
    vals_grid <- expand.grid(comprehenr::to_list(for (var in vnames) c(0,1)))
    colnames(vals_grid) <- vnames
  } else if (type == "2aR50") {
    cat("Testing Example 2a in R50, P_{x}(y|z) is identifiable\n\n")
    adag <- getDAG2aR50()
  } else if (type == "2bR50") {
    cat("Testing Example 2b in R50 (Jiji's example), P_{x}(y|z, w) is identifiable\n\n")
    adag <- getDAG2bR50()
  } else if (type == "2cR50") {
    cat("Testing Example 2c in R50 (backdoor w/ interv), P_{x}(y|z1,z2) is identifiable\n\n")
    adag <- getDAG2cR50()
  } else if (type == "2dR50") {
    cat("Testing Example 2d in R50, P_{x}(y|z1,z2) is identifiable\n\n")
    adag <- getDAG2dR50()
  } else if (type == "2ahytt") {
    cat("Testing Hyttinen's example, P_{x}(y) is identifiable\n\n")
    adag <- getDAG2aHyttinen()
  } else {
    stop("This type has not been implemented yet!")
  }

  if (is.null(vals_grid)) {
    # gets all combinations of adag$x and adag$z with y=1
    vals_grid <- getDefaultValsGrid(adag)
  }
  out <- runIDTest(adag, vals_grid, N)
  return(list(out=out, x=adag$x, y=adag$y, z=adag$z, vals_grid=vals_grid))
}


runIDTest <- function(adag, vals_grid, N=10000, showPlots=FALSE) {
  apag <- getPAG(adag$dagg) # generates the PAG using the perfect oracle
  #plot(apag) # PAG
  amat <- apag@amat
  x <- adag$x
  y <- adag$y
  z <- adag$z

  if (showPlots) {
    dev.new()
    plot(as(t(adag$amat), "graphNEL")) # DAG showing even the latent variables
    dev.new()
    plot(apag) # PAG
    dev.new()
    plotAG(amat) # PAG from the adj matrix
  }

  # done <- FALSE
  # while (!done) {
  #   done <- tryCatch(
  #     {
       # generates a dataset with a random CPT
      dat <- generateDataset(adag$dagg, N, type="binary", alpha_cutoff=0.1, verbose=FALSE)

      #############################
      # Computing P(y | do(x), z) #
      #############################

      adagNEL <- as(t(adag$amat), "graphNEL")
      exprDAG <- causaleffect::causal.effect(y=y, x=x, z=z,
                               G = igraph_from_graphNel(adagNEL, adag$lat),
                               expr = TRUE,
                               simp = TRUE)
      print(paste0("ID from DAG: ", exprDAG))


      retDAG <- causaleffect::causal.effect(y=y, x=x, z=z,
                              G = igraph_from_graphNel(adagNEL, adag$lat),
                              expr = FALSE,
                              simp = TRUE,
                              steps = TRUE)


      other_vnames <- setdiff(colnames(amat), colnames(vals_grid))
      other_vals <- matrix(0, nrow=nrow(vals_grid), ncol=length(other_vnames))
      colnames(other_vals) <- other_vnames
      vals_grid <- cbind(vals_grid, other_vals)


      py.dox.z_DAG <- c()
      for (i in 1:nrow(vals_grid)) {
        curvals <- unlist(vals_grid[i,])
        py.dox.z_DAG <- c(py.dox.z_DAG, estimateInterDistCE(dat, retDAG, vals=curvals))
      }
      py.dox.z_DAG <- cbind(vals_grid, "py_dox.z"=py.dox.z_DAG)


      retPAG <- CIDP(amat, x, y, z)
      retPAG$id
      print(paste0("ID from PAG: ", retPAG$Qexpr[[retPAG$query]]))
      QopElem = retPAG$Qop[[retPAG$query]]
      Qop = retPAG$Qop

      py.dox.z_PAG <- c()
      for (i in 1:nrow(vals_grid)) {
        curvals <- unlist(vals_grid[i,])
        py.dox.z_PAG <- c(py.dox.z_PAG, computeQop(dat, QopElem, Qop, values=curvals)$prob)
      }
      py.dox.z_PAG <- cbind(vals_grid, "py_dox.z"=py.dox.z_PAG)

      differences <- abs(py.dox.z_DAG$py_dox.z - py.dox.z_PAG$py_dox.z)
  #     return(TRUE)
  #     }, error=function(cond) {
  #       message(cond)
  #       return(FALSE)
  #     })
  # }

  # print(py.dox.z_DAG)
  # print(py.dox.z_PAG)
  # print(differences)

  return(list(dag = py.dox.z_DAG, pag = py.dox.z_PAG, diff = differences))
}


runNonIDTests <- function() {
    all_tests <- c()

  #####################
  # Non-ID from PAG 3 #
  #####################

  apag <- getPAG3()
  amat <- apag$amat
  x <- apag$x
  y <- apag$y
  z <- apag$z
  plotAG(amat)

  ret <- CIDP(amat, x, y, z)
  all_tests <- c(all_tests, ret$id)


  #################################
  # Non-ID from PAG Fig.3c in R42 #
  #################################

  apag <- getPAG3cR42()
  amat <- apag$amat
  x <- apag$x
  y <- apag$y
  plotAG(amat)

  ret <- IDP(amat, x, y)
  all_tests <- c(all_tests, ret$id)


  #################################
  # Non-ID from PAG Fig.3a in R42 #
  #################################

  apag <- getPAG3aR42()
  amat <- apag$amat
  x <- apag$x
  y <- apag$y
  plotAG(amat)

  ret <- IDP(amat, x, y)
  all_tests <- c(all_tests, ret$id)

  #################################
  # Non-ID from PAG Fig.3b in R42 #
  #################################

  apag <- getPAG3bR42()
  amat <- apag$amat
  x <- apag$x
  y <- apag$y
  plotAG(amat)

  ret <- IDP(amat, x, y)
  all_tests <- c(all_tests, ret$id)

  #################################
  # Non-ID from PAG Fig.3c in R42 #
  #################################

  apag <- getPAG3cR42()
  amat <- apag$amat
  x <- apag$x
  y <- apag$y
  plotAG(amat)

  ret <- IDP(amat, x, y)
  all_tests <- c(all_tests, ret$id)

  return(list(pass=!any(all_tests), tests=all_tests))
}


dir.create(file.path("./", "out"), showWarnings = FALSE)
nDatasets <- 30
sample_sizes <- c(5*10^3, 10^4, 5*10^4, 10^5, 5*10^5)
tests_diff <- list()
stats <- list()
for (example in c("g1", "g2n", "g2", "2ahytt", "2aR50", "2bR50", "2cR50", "2dR50", "2bR42")) {
  stats[[example]] <- list()
  tests_diff[[example]] <- list()
  for (N in sample_sizes) { # increasing sample sizes
    out_diff <- foreach(nDataset = 1:nDatasets, .combine=rbind) %dopar% {
    cat(paste0("Running example ", example, " with N=", N, "\n\n"))
      #sink(file = paste0("./out/log_", example, "_N", N, "_", nDataset, ".txt"),
      #     append=TRUE)
      test_out <- runTest(type=example, N=N)
      save(test_out, file=paste0("./out/test_out_", example, "_N", N, "_", nDataset, ".RData"))
      test_out$out$diff
    }
    rownames(out_diff) <- 1:nDatasets
    tests_diff[[example]][[as.character(N)]]  <- out_diff

    curmeans <- sapply(as.data.frame(tests_diff[[example]][[as.character(N)]]), mean)
    if (is.null(stats[[example]]$means)) {
      stats[[example]]$means <- curmeans
    } else {
      stats[[example]]$means <- rbind(stats[[example]]$means, curmeans)
    }
    cursds <- sapply(as.data.frame(tests_diff[[example]][[as.character(N)]]), sd)
    if (is.null(stats[[example]]$sds)) {
      stats[[example]]$sds <- cursds
    } else {
      stats[[example]]$sds <- rbind(stats[[example]]$sds, cursds)
    }
    save(tests_diff, file=paste0("./out/tests_diff_", example, "_N", N, ".RData"))
    save(stats, file=paste0("./out/stats_", example, "_N", N, ".RData"))
  }
  rnames <- formatC(sample_sizes, format="f", big.mark=",", digits=0)
  rownames(stats[[example]]$means) <- rnames
  rownames(stats[[example]]$sds) <-rnames

  save(tests_diff, file=paste0("./out/tests_diff_", example, ".RData"))
  save(stats, file=paste0("./out/stats_", example, ".RData"))
}
save(stats, file=paste0("./out/stats.RData"))
plotResults = FALSE

####################
# Plotting Results #
####################
library(ggplot2)

getVLabels <- function(vnames, x, y, z, vals_grid) {
  vlabels <- c()
  for (v in vnames) {
    xvalues <- paste(paste(x, "=", vals_grid[v, x], sep=""), collapse=",")
    yvalues <- paste(paste(y, "=", vals_grid[v, y], sep=""), collapse=",")
    if (!is.null(z)) {
      zvalues <- paste(paste(z, "=", vals_grid[v, z], sep=""), collapse=",")
      vlabel <- paste("P(", yvalues, "| do(", xvalues, "), ", zvalues, ")", sep="")
    } else {
      vlabel <- paste("P(", yvalues, "| do(", xvalues, "))", sep="")
    }

    vlabels <- c(vlabels, vlabel)
  }
  names(vlabels) = 1:length(vnames)
  return(vlabels)
}

getSimSummary <- function(curmeans, cursds) {
  sim_summary <- data.frame()
  for (v in 1:ncol(cursds)) {
    cur_v_sim <- c()
    for (r in 1:nrow(cursds)) {
      ymin = curmeans[r,v] - cursds[r,v]
      ymax = curmeans[r,v] + cursds[r,v]
      ymean = curmeans[r,v]
      cur_v_sim <-  rbind(cur_v_sim, cbind(ymin, ymax, ymean))
    }
    nextrows <- (nrow(sim_summary)+1):(nrow(sim_summary)+nrow(cursds))
    sim_summary[nextrows,c("ymin", "ymax", "ymean")] <- cur_v_sim
    sim_summary[nextrows,"N"] <- factor(rnames, levels=rnames)
    sim_summary[nextrows,"V"] <- v
  }
  sim_summary[,"V"] <- as.factor(sim_summary[,"V"])
  return(sim_summary)
}

plotErrorBars <- function(sim_summary, vlabels, saveFile, nrow=4, width=9, height=5) {
  aplot <- ggplot(sim_summary, aes(x = N, y = ymean)) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = ymin, ymax = ymax)) +
    ggtitle("Average absolute difference in causal effect given causal diagram and PAG") +
    theme(plot.title = element_text(size = 16, face = "bold"),
          axis.text.x = element_text(size=10, angle = 30),
          axis.text.y = element_text(size = 10),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          strip.text = element_text(size = 10)) +
    xlab("Number of observations in each dataset") +
    ylab("Average absolute difference from 30 datasets") +
    facet_wrap( ~ V, nrow=nrow, labeller = labeller(V = vlabels ))
  # labs(title="hwy vs displ", caption = "Source: mpg", subtitle="Ggplot2 - Faceting - Multiple plots in one figure")  # Shared scales

  if (!is.null(saveFile)) {
    pdf(saveFile, width = width, height = height)
    print(aplot)
    dev.off()
  } else{
    return(aplot)
  }
}

sample_sizes <- c(5*10^3, 10^4, 5*10^4, 10^5, 5*10^5)
ex_folder <- "./out/"
if (plotResults) {
  for (example in c("g1", "g2", "2bR42", "2ahytt", "2aR50", "2bR50", "2cR50", "2dR50")) {
    load(file=paste0(ex_folder, "test_out_", example, "_N", sample_sizes[1], "_1.RData"))
    vals_grid <- test_out$vals_grid
    x <- test_out$x
    y <- test_out$y
    z <- test_out$z

    load(file=paste0(ex_folder, "stats_", example, ".RData"))
    load(file=paste0(ex_folder, "tests_diff_", example, ".RData"))

    rnames <- formatC(sample_sizes, format="f", big.mark=",", digits=0)
    rownames(stats[[example]]$means) <- rnames
    rownames(stats[[example]]$sds) <-rnames

    cursds <- stats[[example]]$sds
    curmeans <- stats[[example]]$means
    vnames <- 1:ncol(cursds)

    vlabels <- getVLabels(vnames, x, y, z, vals_grid)
    sim_summary <- getSimSummary(curmeans, cursds)

    nrow <- max(2,ceiling(sqrt(length(vlabels))))
    ncol <- ceiling(length(vlabels)/nrow)
    plotErrorBars(sim_summary, vlabels,
                  saveFile = paste0("./plot_", example, ".pdf"), nrow=nrow, height=max(nrow*1.5, 3.5), width=ncol*2.5)
  }
}
