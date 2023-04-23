# library(methods)
library(comprehenr)

getPAG <- function(dag, verbose = FALSE) {
  valR <- FALSE
  while(!valR) {
    R <- dagitty::impliedCovarianceMatrix(dag, b.default = NULL, b.lower = -0.6,
                                          b.upper = 0.6, eps = 1, standardized = TRUE)

    R <- round(R, 14)
    valR <- matrixcalc::is.symmetric.matrix(R) &&
      matrixcalc::is.positive.definite(R, tol=1e-8)
    if(verbose)
      cat(paste0("valR=", valR, "\n"))
  }

  latR <- R
  suffStat = list(C = latR, n = 10^9)

  true.pag <- pcalg::fci(suffStat,
                         indepTest = pcalg::gaussCItest, #p = ncol(-latR),
                         labels= colnames(suffStat$C), alpha = 0.9999)
  #plot(true.pag)
  return(true.pag)
}

#colSds <- function (x, na.rm=FALSE) apply(X=x, MARGIN=2, FUN=sd, na.rm=na.rm)

# Generates Obs. dataset according to the DAG #
# adag is either a dagitty or a bn.fit object
# TODO: varlevels can be a named array with the levels of each variable in adag
generateDataset <- function(adag, N, type="binary", alpha_cutoff = 0.2, varlevels=3, checkLocTests=TRUE, forcePosit=TRUE, verbose=FALSE) {
  # TODO: local tests with conditional dependencies!

  if (!(any(class(adag) == "bn.fit") || any(class(adag) == "dagitty")))  {
    stop("adag must be a dagitty or a bn.fit object")
  }

  setobs <- FALSE
  while (!setobs) {
    setobs <- tryCatch(
      {
        if (any(class(adag) == "bn.fit")) {
          obs.dat <- bnlearn::rbn(adag, N)
          obs.dat <- as.data.frame(sapply(obs.dat, function(col) as.numeric(col)-1))
          amat <- t(bnlearn::amat(adag)) # it needs to be transposed to meet pcalg notation
          #plot(as(t(amat), "graphNEL")) # DAG
          adagg <-  pcalg::pcalg2dagitty(amat, colnames(amat),type="dag")
          if (type =="binary" || type == "categorical") {
            lt <- dagitty::localTests(adagg, obs.dat, type="cis.chisq")
            #print(lt)
          }
        }
        else if (any(class(adag) == "dagitty")) {
          dagg <- adag
          if(type == "binary") {
            obs.dat <- dagitty::simulateLogistic(adag, N=N, verbose=FALSE)
            obs.dat <- as.data.frame(sapply(obs.dat, function(col) as.numeric(col)-1))
            lt <- dagitty::localTests(adag, obs.dat, type="cis.chisq")
            min(lt$p.value)
            #print(lt)
          } else if (type == "linear") {
            obs.dat <- dagitty::simulateSEM(adag, N=N)
            lt <- dagitty::localTests(adag, obs.dat, type="cis")
          } else if (type == "categorical") {
            obs.dat <- simulateSEM(adag, N=N)
            obs.dat <- as.data.frame(lapply(obs.dat, function(col) cut(col, varlevels)))
            lt <- dagitty::localTests(adag, obs.dat, type="cis.chisq")
            #print(lt)
          }
        }

        alpha <- min(lt$p.value)
        if (verbose)
          print(alpha)

        ispos <- TRUE
        if (checkLocTests && alpha < alpha_cutoff) {
          cat("unfaithful dataset...\n")
          ispos <- FALSE
        }
        if (ispos && forcePosit) {
          vnames <- colnames(obs.dat)
          vals_grid <- expand.grid(comprehenr::to_list(for (var in vnames) c(0,1)))
          colnames(vals_grid) <- vnames
          i = 1
          while (ispos == TRUE && i <= nrow(vals_grid)) {
            curcond <- getCondProb(obs.dat, var=vnames, cond=NULL, vals=vals_grid[i,])
            if (curcond == 0) {
              cat("non-positive dataset...\n")
              cat(paste0(vnames, "=", vals_grid[i,], collapse=";"), "\n")
              ispos = FALSE
            }
            i = i+1
          }
        }
        ispos
    }, error=function(cond) {
      message(cond)
      return(FALSE)
    })
  }

  #head(obs.dat)

  # if (binary_data) {
  #   Xname <- exposures(dagg)
  #   Yname <- outcomes(dagg)
  #
  #   # P(y=1|X=0)
  #   py1x0 <- length(intersect( which(obs.dat[,Yname] == 1) , which(obs.dat[,Xname] == 0) ))/N
  #   px0 <- length(which(obs.dat[,Xname] == 0))/N
  #   py1.x0 <- py1x0 / px0
  #
  #   # P(y=1|X=1)
  #   py1x1 <- length(intersect( which(obs.dat[,Yname] == 1) , which(obs.dat[,Xname] == 1) ))/N
  #   px1 <- length(which(obs.dat[,Xname] == 1))/N
  #   py1.x1 <- py1x1 / px1
  # }

  return(obs.dat)
}


getCondProb <- function(dat, var, cond, vals) {
  if (all(sapply(dat, function(col) is.character(col)))) {
    values <- paste0("\"",  vals[c(var, cond)], "\"")
    values_den <- paste0("\"",  vals[cond], "\"")
  } else if (all(sapply(dat, function(col) is.numeric(col)))) {
    values <- vals[c(var, cond)]
    values_den <- vals[cond]
  } else {
    stop("The type of all variables in dat should be either numeric or character")
  }

  # print(paste0("var: ", paste0(var, collapse=",")))
  # print(paste0("cond:  ", paste0(cond, collapse=",")))
  # print(paste0("vals: ", paste0(names(vals), collapse=","), " ={", paste0(vals, collapse=","), "}"))

  condition_num <- paste(c(var, cond), "==", values, collapse=" & ")
  # print(paste0("condition_num: ", condition_num))

  num <- dim(subset(dat, subset=eval(parse(text=condition_num))))[1]

  if (length(cond) > 0) {
    condition_den <-  paste(paste(c(cond), "==", values_den, collapse=" & "))
    #print(paste0("condition_den: ", condition_den))
    den <- dim(subset(dat, subset=eval(parse(text=condition_den))))[1]
  } else {
    den <- dim(dat)[1]
  }

  if (den <= 0){
    den <- 1
  }

  # print(paste0("num prob: ", num))
  # print(paste0("den prob: ", den))
  # print(paste0("frac prob: ", num/den))

  prob <-  num /den
  return(prob)
}



computeQopNone <- function(dat, QopElem, Qop, values, jointProbTable) {
  if (length(QopElem$param$interv) > 0) {
    joint_out <- computeQop(dat,
                            Qop[[paste0(QopElem$param$interv, collapse=",")]],
                            Qop, values=values, jointProbTable)
    jointProbTable <- joint_out$jpt
    cur_joint <- joint_out$prob
  } else {
    jpt_cond <- paste0(paste0(names(values), " == ", values), collapse=" & ")
    cur_jpt <- subset(jointProbTable, subset=eval(parse(text=jpt_cond)))

    if (dim(cur_jpt)[1] > 0) {
      cur_joint <- cur_jpt[1,"prob"]
    } else {
      cur_joint <- getCondProb(dat, var=names(values), cond=NULL, vals=values)
      tab_valnames <- colnames(jointProbTable)
      tab_valnames <- tab_valnames[-length(tab_valnames)]
      jointProbTable[nrow(jointProbTable)+1,] <- c(values[tab_valnames], "prob"=cur_joint)
    }
  }
  return(list(prob=cur_joint, jpt=jointProbTable))
}

getLevels <- function(dat, varname) {
  var <- dat[,varname]

  if (is.factor(var)) {
    return(levels(var))
  } else if (is.numeric(var)) {
    return(sort(unique(var)))
  } else {
    stop(paste0("Type of ", varname, ", ", class(var), ", is not supported yet"))
  }
}


computeOpFracCond <- function(dat, QopElem, Qop, values, jointProbTable) {
  # prob = param$prob / \sum_den.sumset param$prob
  if (length(QopElem$param$prob) > 0) { # vs obs #TODO check
    num_out <- computeQop(dat, Qop[[QopElem$param$prob]], Qop, values, jointProbTable)
    jointProbTable <- num_out$jpt
    num <- num_out$prob
  } else {
    jpt_cond <- paste0(paste0(names(values), " == ", values), collapse=" & ")
    cur_jpt <- subset(jointProbTable, subset=eval(parse(text=jpt_cond)))

    if (dim(cur_jpt)[1] > 0) {
      num <- cur_jpt[1,"prob"]
    } else {
      num <- getCondProb(dat, var=names(values), cond=NULL, vals=values)
      tab_valnames <- colnames(jointProbTable)
      tab_valnames <- tab_valnames[-length(tab_valnames)]
      jointProbTable[nrow(jointProbTable)+1,] <- c(values[tab_valnames], "prob"=num)
    }
  }

  sumset <- QopElem$param$den.sumset
  sumset_values <- expand.grid(comprehenr::to_list(for (var in sumset) getLevels(dat, var)))
  colnames(sumset_values) <- sumset

  den <- 0
  den_values <- values
  for (irow in 1:nrow(sumset_values)) {
    den_values[sumset] <- sumset_values[irow,]
    if (length(QopElem$param$prob) > 0) {
      den_prob_out <- computeQop(dat, Qop[[QopElem$param$prob]], Qop, values=den_values, jointProbTable)
      jointProbTable <- den_prob_out$jpt
      den_prob <- den_prob_out$prob
    } else {
      jpt_cond <- paste0(paste0(names(den_values), " == ", den_values), collapse=" & ")
      cur_jpt <- subset(jointProbTable, subset=eval(parse(text=jpt_cond)))

      if (dim(cur_jpt)[1] > 0) {
        den_prob <- cur_jpt[1,"prob"]
      } else {
        den_prob <- getCondProb(dat, var=names(den_values), cond=NULL, vals=den_values)
        tab_valnames <- colnames(jointProbTable)
        tab_valnames <- tab_valnames[-length(tab_valnames)]
        jointProbTable[nrow(jointProbTable)+1,] <- c(den_values[tab_valnames], "prob"=den_prob)
      }
    }
    den <- den + den_prob
  }

  if (num == 0 && den == 0) {
    print("both numerator and denominator are zero!")
    final_prob <- 0
  } else if (den == 0) {
    stop("denominator is zero while numerator is not!")
  } else {
    final_prob <- num/den
  }

  return(list(prob=final_prob, jpt=jointProbTable))
}

computeOpSumset <- function(dat, QopElem, Qop, values, jointProbTable) {
  sumset <- QopElem$param$sumset
  sumset_values <- expand.grid(comprehenr::to_list(for (var in sumset) getLevels(dat, var)))
  colnames(sumset_values) <- sumset

  joint_prob <- 0
  joint_values <- values
  for (irow in 1:nrow(sumset_values)) {
    joint_values[sumset] <- sumset_values[irow,]
    joint_values <- unlist(joint_values)

    if (length(QopElem$param$interv) > 0) {
      cur_joint_out <- computeQop(dat,
                                  Qop[[paste0(QopElem$param$interv, collapse=",")]],
                                  Qop, values=joint_values, jointProbTable)
      jointProbTable <- cur_joint_out$jpt
      cur_joint <- cur_joint_out$prob
    } else {
      jpt_cond <- paste0(paste0(names(joint_values), " == ", joint_values), collapse=" & ")
      cur_jpt <- subset(jointProbTable, subset=eval(parse(text=jpt_cond)))
      if (dim(cur_jpt)[1] > 0) {
        cur_joint <- cur_jpt[1,"prob"]
      } else {
        cur_joint <- getCondProb(dat, var=names(joint_values), cond=NULL, vals=joint_values)
        tab_valnames <- colnames(jointProbTable)
        tab_valnames <- tab_valnames[-length(tab_valnames)]
        jointProbTable[nrow(jointProbTable)+1,] <- c(joint_values[tab_valnames], "prob"=cur_joint)
      }
    }
    joint_prob <- joint_prob + cur_joint
  }

  return(list(prob=joint_prob, jpt=jointProbTable))
}


computeOpProp7 <- function(dat, QopElem, Qop, values, jointProbTable) {
  num.prod1_out = computeQop(dat,
                             Qop[[paste0(QopElem$param$num.prod1, collapse=",")]],
                             Qop, values=values, jointProbTable)
  jointProbTable <- num.prod1_out$jpt
  num.prod1 <- num.prod1_out$prob

  num.prod2_out = computeQop(dat,
                             Qop[[paste0(QopElem$param$num.prod2, collapse=",")]],
                             Qop, values=values, jointProbTable)
  jointProbTable <- num.prod2_out$jpt
  num.prod2 <- num.prod2_out$prob

  den_out = computeQop(dat,
                       Qop[[paste0(QopElem$param$den, collapse=",")]],
                       Qop, values=values, jointProbTable)
  jointProbTable <- den_out$jpt
  den <- den_out$prob

  if (den == 0) {
    stop("denominator is equal to zero!")
  }

  return(list(prob=((num.prod1 * num.prod2) / den), jpt=jointProbTable))
}


computeOpProp6 <- function(dat, QopElem, Qop, values, jointProbTable) {
  # prob =  param$interv(num.var) / param$interv(den.var | den.cond)

  all_vars <- c(QopElem$param$interv, QopElem$param$num.var)
  all_values <- values[all_vars]
  if (length(all_values) != length(values)) {
    # Does all_values always equal to values?
    stop("Check all_values!")
  }

  #######################
  # computing numerator #
  #######################
  if (length(QopElem$param$interv) > 0) {
    num_joint_prob_out <- computeQop(dat, Qop[[paste0(QopElem$param$interv, collapse=",")]], Qop, values=all_values, jointProbTable)
    jointProbTable <- num_joint_prob_out$jpt
    num_joint_prob <- num_joint_prob_out$prob
  } else {
    jpt_cond <- paste0(paste0(names(all_values), " == ", all_values), collapse=" & ")
    cur_jpt <- subset(jointProbTable, subset=eval(parse(text=jpt_cond)))

    if (dim(cur_jpt)[1] > 0) {
      num_joint_prob <- cur_jpt[1,"prob"]
    } else {
      num_joint_prob <- getCondProb(dat, var=names(all_values), cond=NULL, vals=all_values)
      tab_valnames <- colnames(jointProbTable)
      tab_valnames <- tab_valnames[-length(tab_valnames)]
      jointProbTable[nrow(jointProbTable)+1,] <- c(all_values[tab_valnames], "prob"=num_joint_prob)
    }
  }

  #########################
  # computing denominator #
  #########################

  #########################################
  # Joint distribution in the denominator #
  #########################################
  den_joint_vars <- c(QopElem$param$den.var, QopElem$param$den.cond)
  if (all(QopElem$param$num.var %in% den_joint_vars)) {
    den_joint_prob <- num_joint_prob
  } else {
    # getting the joint distributions of the denominator by
    # marginalizing out the variables that only appear in the numerator
    sumset_den <- setdiff(QopElem$param$num.var, den_joint_vars)
    sumset_den_values <- expand.grid(comprehenr::to_list(for (var in sumset_den) getLevels(dat, var)))
    colnames(sumset_den_values) <- sumset_den

    den_joint_prob <- 0
    den_joint_values <- all_values
    for (irow in 1:nrow(sumset_den_values)) {
      den_joint_values[sumset_den] <- sumset_den_values[irow,]
      if (all(unlist(den_joint_values) == all_values)) {
        cur_den_joint <- num_joint_prob
      } else {
        if (length(QopElem$param$interv) > 0) {
          cur_den_joint_out <- computeQop(dat, Qop[[paste0(QopElem$param$interv, collapse=",")]], Qop, values=den_joint_values, jointProbTable)
          jointProbTable <- cur_den_joint_out$jpt
          cur_den_joint <- cur_den_joint_out$prob
        } else {
          jpt_cond <- paste0(paste0(names(den_joint_values), " == ", den_joint_values), collapse=" & ")
          cur_jpt <- subset(jointProbTable, subset=eval(parse(text=jpt_cond)))

          if (dim(cur_jpt)[1] > 0) {
            cur_den_joint <- cur_jpt[1,"prob"]
          } else {
            cur_den_joint <- getCondProb(dat, var=names(den_joint_values), cond=NULL, vals=den_joint_values)
            tab_valnames <- colnames(jointProbTable)
            tab_valnames <- tab_valnames[-length(tab_valnames)]
            jointProbTable[nrow(jointProbTable)+1,] <- c(den_joint_values[tab_valnames], "prob"=cur_den_joint)
          }
        }
      }
      den_joint_prob <- den_joint_prob + cur_den_joint
    }
  }

  ###################################################################
  #  marginal distribution of the variables in the conditioning set #
  ###################################################################
  den_cond_vars <- QopElem$param$den.cond
  sumset_cond <- setdiff(QopElem$param$num.var, den_cond_vars)
  sumset_cond_values <- expand.grid(comprehenr::to_list(for (var in sumset_cond) getLevels(dat, var)))
  colnames(sumset_cond_values) <- sumset_cond

  den_cond_prob <- 0
  den_cond_values <- all_values
  for (irow in 1:nrow(sumset_cond_values)) {
    den_cond_values[sumset_cond] <- sumset_cond_values[irow,]
    if (all(unlist(den_cond_values) == all_values)) {
      cur_den_cond <- num_joint_prob
    } else {
      if (length(QopElem$param$interv) > 0) {
        cur_den_cond_out <- computeQop(dat, Qop[[paste0(QopElem$param$interv, collapse=",")]], Qop, values=den_cond_values, jointProbTable)
        jointProbTable <- cur_den_cond_out$jpt
        cur_den_cond <- cur_den_cond_out$prob
      } else {
        jpt_cond <- paste0(paste0(names(den_cond_values), " == ", den_cond_values), collapse=" & ")
        cur_jpt <- subset(jointProbTable, subset=eval(parse(text=jpt_cond)))

        if (dim(cur_jpt)[1] > 0) {
          cur_den_cond <- cur_jpt[1,"prob"]
        } else {
          cur_den_cond <- getCondProb(dat, var=names(den_cond_values), cond=NULL, vals=den_cond_values)
          tab_valnames <- colnames(jointProbTable)
          tab_valnames <- tab_valnames[-length(tab_valnames)]
          jointProbTable[nrow(jointProbTable)+1,] <- c(den_cond_values[tab_valnames], "prob"=cur_den_cond)
        }
      }
    }
    den_cond_prob <- den_cond_prob + cur_den_cond
  }

  if (den_joint_prob == 0 && den_cond_prob == 0) {
    print("!!!!!!!both den_joint_prob and den_cond_prob are 0!!!!!!")
    den_prob <- 0
  } else if (den_cond_prob == 0) {
    stop("ERROR: den_cond_prob is zero but den_joint_prob is not!!!!")
  } else {
    den_prob <- den_joint_prob/den_cond_prob
  }

  if (num_joint_prob == 0 && den_prob == 0) {
    print("!!!!!!!both num_joint_prob and den_prob are 0!!!!!!")
    final_prob <- 0
  } else if (den_prob == 0) {
    stop("ERROR: den_prob is zero but num_joint_prob is not!!!!")
  } else {
    final_prob <- num_joint_prob/den_prob
  }

  return(list(prob=final_prob, jpt=jointProbTable))
}

# jointProbTable is a table with the computed joint probabilities
computeQop <- function(dat, QopElem, Qop, values, jointProbTable = NULL) {
  #print(paste0("computing operation ", QopElem$type))

  if (is.null(jointProbTable)) {
    jointProbTable <- data.frame(matrix(ncol=length(values)+1, nrow = 0))
    colnames(jointProbTable) <- c(names(values), "prob")
  }

  if (QopElem$type == "none") {
    return(computeQopNone(dat, QopElem, Qop, values, jointProbTable))
  } else if (QopElem$type == "frac_cond") {
    return(computeOpFracCond(dat, QopElem, Qop, values, jointProbTable))
  } else if (QopElem$type == "sumset") {
    return(computeOpSumset(dat, QopElem, Qop, values, jointProbTable))
  } else if (QopElem$type == "prop7") {
    return(computeOpProp7(dat, QopElem, Qop, values, jointProbTable))
  } else if (QopElem$type == "prop6") {
    return(computeOpProp6(dat, QopElem, Qop, values, jointProbTable))
  }
}
