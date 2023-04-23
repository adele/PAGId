rm(list=ls())

library(pcalg)
library(dagitty)
library(igraph)
library(causaleffect)
library(R.utils)
library(PAGId)
source("PAGInstances.R")
source("ProbUtils.R")
source("EstimateCausalEffect.R")

source('../Solver PAG-ID/docalcwograph/R/load.R')

loud()

#################
# Aux Functions #
#################


# Converts G, GE to pcalg amat
amatFromGGe <- function(G, Ge, vnames) {
  amat <- G
  colnames(amat) <- rownames(amat) <- vnames
  #plot(as(t(amat), "graphNEL"))
  latv <- c()
  for (i in 1:nrow(Ge)) {
    confv <- which(Ge[,i] == 1 & Ge[i,] == 1)
    if (length(confv) > 0) {
      for (j in 1:length(confv)) {
        ivar <- vnames[i]
        jvar <- vnames[confv[j]]
        curlatv <- paste0("u", ivar, jvar)
        curlatv2 <- paste0("u", jvar, ivar)
        curnames <- colnames(amat)
        if (any(c(curlatv, curlatv2) %in% curnames))
          next
        latv <- c(latv, curlatv)
        amat <- cbind(amat, 0)
        amat <- rbind(amat, 0)
        colnames(amat) <- rownames(amat) <- c(curnames, curlatv)

        amat[curlatv, ivar] <- 0; amat[ivar, curlatv] <- 1; # curlatv -> i
        amat[curlatv, jvar] <- 0; amat[jvar, curlatv] <- 1; # curlatv -> j
      }
    }
  }
  #plot(as(t(amat), "graphNEL"))
  return(list(amat=amat, lat=latv))
}

rundocalcwograph <- function(L, noind=FALSE, timeout=300) {
  ret_docalcwog <-
    tryCatch(
    {
      formular <- R.utils::withTimeout(
        {
          #if (noind) {
          #  setdocalc.noind(L,fix_order=FALSE)
          #} else {
            setdocalc(L,fix_order=FALSE)
          #}
        },
        onTimeout="error",
        timeout=timeout,
        cpu=timeout,
        elapsed=timeout)

      if (!is.null(formular) && !any(is.na(formular$formulas))) {
        print(paste0("ID from PAG using docalcwograph"))
        print(formular$formulas)
        ret <- TRUE
      } else {
        print(paste0("Non-ID from PAG using docalcwograph"))
        ret <- FALSE
      }
      list(id=ret, formular=formular)
    }, TimeoutException = function(ex) {
      message("Timeout in docalcwograph")
      list(id=FALSE, formular=list(times="timeout"))
    }, error=function(err) {
      message(paste0("Cathing error in docalcwograph PAG ID: ", err))
      list(id=FALSE, formular=NULL)
    })

  return(ret_docalcwog)
}


runCIDPgivenPAG <- function(amat_pag, x, y, z) {
  times_cidp <- 0-proc.time()[3]
  retPAG <- CIDP(amat_pag, x, y, z, verbose=FALSE)
  times_cidp <- times_cidp + proc.time()[3]

  if (retPAG$id) {
    id_CIDP <- TRUE
    print(paste0("ID from PAG using CIDP: ", retPAG$Qexpr[[retPAG$query]]))
    #retPAG$Qexpr
  } else {
    id_CIDP <- FALSE
    print(paste0("Non-ID from PAG using CIDP"))
  }
  return(list(id=id_CIDP, time_cidp=times_cidp, ret=retPAG))
}


runPAGandCIDP <- function(amat, x, y, z) {
  lati <- which(substr(colnames(amat), 1, 1) == "u")
  if (length(lati) > 0) {
    lat <- colnames(amat)[lati]
  } else {
    lat <- c()
  }
  adag <- pcalg::pcalg2dagitty(amat, colnames(amat),type="dag")
  dagitty::latents(adag) <- lat
  dagitty::exposures(adag) <- x
  dagitty::outcomes(adag) <- y
  #plot(adag)

  ###########################################
  # Checking if the effect if identifiable  #
  # from the PAG using CIDP                 #
  ###########################################

  apag <- NULL
  pag_timeout <- 600 #5 minutes
  times_CIDP <- list()
  ret_apag <-
    tryCatch(
      {
        apag <- R.utils::withTimeout(
          {
            times_CIDP$pag <- 0-proc.time()[3]
            apag_ret <- getPAG(adag)
            times_CIDP$pag <- times_CIDP$pag+proc.time()[3] # time in seconds
            apag_ret
          },
          onTimeout="error",
          timeout=pag_timeout,
          cpu=pag_timeout,
          elapsed=pag_timeout)
        list(err_timeout=FALSE, apag=apag)
      }, TimeoutException = function(ex) {
        message("Timeout in generating a PAG")
        list(err_timeout=TRUE, apag=NULL)
      }, error=function(err) {
        message(paste0("Unexpected error while generating a PAG: ", err))
        list(err_timeout=FALSE, apag=NULL)
      })

  apag <- ret_apag$apag
  if (is.null(apag))
    return(NULL)

  # if (plotg) {
  #   graphics.off()
  #   plot(apag)
  # }

  amat_pag <- apag@amat
  times_CIDP$cipd <- 0-proc.time()[3]
  retPAG <- CIDP(amat_pag, x, y, z, verbose=TRUE)
  times_CIDP$cipd <- times_CIDP$cipd+proc.time()[3]

  if (retPAG$id) {
    id_CIDP <- TRUE
    print(paste0("ID from PAG using CIDP: ", retPAG$Qexpr[[retPAG$query]]))
    #retPAG$Qexpr
  } else {
    id_CIDP <- FALSE
    print(paste0("Non-ID from PAG using CIDP"))
  }
  return(list(id=id_CIDP, times=times_CIDP, amat_pag=amat_pag, ret=retPAG))
}

###################################
# Generating docalcwograph object #
###################################

generateRandomGraph <- function(n, nedges_perc) {
  L <- list()

  nedges <- nedges_perc * (2*(n*(n-1)/2))

  L$n<-n
  cat("Generating random graph:\n\n\n")
  M<-random.do.graph(n,nedge=nedges,nconf=NA,ancforce=TRUE,model_index=NA )

  G<-M$G
  Ge<-M$Ge

  L$U<-2 # Y
  L$J<-1 # X
  L$C<-c(1,3) # X, Z
  vnames <-  c("X", "Y", paste0("V", 1:(L$n-2)))
  L$vnames <- vnames

  #cat("Generating discrete random model:\n\n\n")
  L$M<-randomModel.discrete(G,Ge)

  uvec<-rep(0,L$n)
  uvec[L$U]<-1
  cvec<-rep(0,L$n)
  cvec[L$C]<-1
  jvec<-rep(0,L$n)
  jvec[L$J]<-1

  L$uset<-bin.to.dec(rev(uvec))
  L$cset<-bin.to.dec(rev(cvec))
  L$jset<-bin.to.dec(rev(jvec))

  L$start = paste('p(',paste(L$U,collapse=','),'|',paste(L$C,collapse=','),'||',paste(L$J,collapse=','),')',sep='')
  L$query =  paste('p(',L$uset,',',L$cset,',',L$jset,')',sep='')
  return(L)
}

########################################################
# Generating Random DAGs using Hyttinen et al.'s code. #
########################################################

plotg = FALSE
nedges_perc = 0.4
idcase = TRUE
runCIDP = TRUE
rundocalcwog = FALSE

for (n in 13:15) {
  i = 1
  while(i <= 30) {
    if (n == 13 && i <= 6) {
     i = i+1
     next
    }
    L <- generateRandomGraph(n, nedges_perc)

    ###########################################
    # Checking if the effect if identifiable  #
    # from the DAG using causaleffect package #
    ###########################################

    cat("Cheking ID from DAG:\n\n\n")

    ig <- m.to.igraph(L$M, L$vnames)
    if (plotg) {
      graphics.off()
      plot(ig)
    }

    id_dag <- tryCatch(
      {
        exprDAG <- causaleffect::causal.effect(y=L$vnames[L$U],
                                               x=L$vnames[L$J],
                                               z=L$vnames[setdiff(L$C, L$J)],
                                               G = ig,
                                               expr = TRUE,
                                               simp = TRUE)

        retDAG <- NULL
        retDAG <- causaleffect::causal.effect(y=L$vnames[L$U],
                                               x=L$vnames[L$J],
                                               z=L$vnames[setdiff(L$C, L$J)],
                                               G = ig,
                                               expr = FALSE,
                                               simp = TRUE, steps = TRUE)

        ret <- FALSE
        if (!is.null(retDAG) && retDAG$id) { # id
          ret <- TRUE
          if (idcase && !any(c(retDAG$P$sum,retDAG$P$product,retDAG$P$fraction))) {
            print(paste0("ID from DAG: ", exprDAG))
            ret <- FALSE
          }
        }
        ret
      }, error=function(err) {
        message(paste0("Cathing error in DAG ID: ", err))
        FALSE
      })

    if (id_dag == idcase) {
      cat(paste0("ID from DAG:", id_dag, "\n\n\n"))
      amat_out <- amatFromGGe(L$M$G, L$M$Ge, L$vnames)
      amat <- amat_out$amat

      if (plotg) {
        graphics.off()
        agraph <- as(t(amat), "graphNEL")
        plot(agraph)
      }

      ##################################
      # Getting the PAG e running CIDP #
      ##################################
      x <- L$vnames[L$J]
      y <- L$vnames[L$U]
      z <- L$vnames[setdiff(L$C, L$J)]

      if (runCIDP) {
        retCIDP <- runPAGandCIDP(amat, x, y, z) #runTruePAGandCIDP
        id_CIDP <- retCIDP$id
        times_CIDP <- retCIDP$times
        amat_pag <- retCIDP$amat_pag
      } else {
        amat_pag <- NULL
        times_CIDP <- NULL
        id_CIDP <- "PENDING"
      }

      if (idcase && !id_CIDP)
        next

      ###########################################
      # Checking if the effect if identifiable  #
      # from the PAG using CIDP                 #
      ###########################################

      if (rundocalcwog) {
        ret_docalcwog <- rundocalcwograph(L, noind = !idcase)
        id_docalcwog <- ret_docalcwog$id
        formular <- ret_docalcwog$formular
      } else {
        formular <- NULL
        id_docalcwog <- "PENDING"
      }

      print(paste0("Saving objects... CIDP_", id_CIDP, "_docalcwog_", id_docalcwog, "_", i))
      saveobj = list(L=L, ig=ig, amat_dag=amat, amat_pag=amat_pag, x=x, y=y, z=z, times_CIDP=times_CIDP, formular=formular)
      save(saveobj, file=paste0("./random/ins_", n, "_CIDP_", id_CIDP, "_docalcwog_", id_docalcwog, "_", i, ".RData"))
      i = i+1
    }
  }
}
`
dirname <- "./random_nonid/"
rundocalcwog <- TRUE
runCIDP <- FALSE
for (n in 5:20) {
  for (i in 1:30) {
    if (n == 5 && i <= 26) {
      next
    }
    saveobj <- NULL
    formular <- NULL

    filename <- paste0(dirname, "ins_", n, "_CIDP_FALSE_docalcwog_FALSE_", i, ".RData")
    if (!file.exists((filename))) {
      filename <- paste0(dirname, "ins_", n, "_CIDP_FALSE_docalcwog_PENDING_", i, ".RData")
    }

    if (file.exists((filename))) {
      load(filename)

      amat <- saveobj$amat_dag
      x <- saveobj$x
      y <- saveobj$y
      z <- saveobj$z

      if (runCIDP) {
        retCIDP <- runPAGandCIDP(amat, x, y, z)
        id_CIPD <- retCIDP$id
        times_CIDP <- retCIDP$times
        saveobj$times_CIDP <- times_CIDP
        id_CIDP <- retCIDP$id
      } else {
        id_CIDP <- FALSE
      }

      if (rundocalcwog) {
        L <- saveobj$L
        ret_docalcwog <- rundocalcwograph(L, timeout = 1200) #10 minutes of timeout
        id_docalcwog <- ret_docalcwog$id
        formular <- ret_docalcwog$formular
        if (!is.null(formular) && length(formular$times) == 1)
          id_docalcwog <- "TIMEOUT"
        saveobj$formular <- formular
      } else {
        id_docalcwog <- "PENDING"
      }
      save(saveobj, file=paste0("./random/ins_", n, "_CIDP_", id_CIDP, "_docalcwog_", id_docalcwog, "_", i, ".RData"))
    } else {
      cat(paste0("ERROR - file not found: ", filename))
    }
  }
}



# 1143 10   3 CIDP ID   17.250
# 1173 11   3 CIDP ID  170.951
# 1179 11   9 CIDP ID 1258.776
# 1181 11  11 CIDP ID   84.062
# 1205 12   5 CIDP ID   75.054
# 1206 12   6 CIDP ID   44.462
# 1227 12  27 CIDP ID   55.573
n = 10
i = 3

idcase <- FALSE
if (idcase) {
  curdir <- "./random_id/"
} else {
  curdir <- "./random_nonid/"
}
for (n in 5:12) {
  for (i in 1:30) {
    filename <- paste0(curdir, "ins_", n, "_CIDP_", idcase, "_docalcwog_TRUE_", i, ".RData")
    id_docalcwog <- TRUE
    id_CIDP <- idcase
    if (!file.exists((filename))) {
      filename <- paste0(curdir, "ins_", n, "_CIDP_", idcase, "_docalcwog_FALSE_", i, ".RData")
      id_docalcwog <- FALSE
      id_CIDP <- idcase
    }

    if (!file.exists((filename))) {
      filename <- paste0(curdir, "ins_", n, "_CIDP_", idcase, "_docalcwog_PENDING_", i, ".RData")
      id_docalcwog <- "PENDING"
      id_CIDP <- idcase
    }

    if (!file.exists((filename))) {
      filename <- paste0(curdir, "ins_", n, "_CIDP_", idcase, "_docalcwog_TIMEOUT_", i, ".RData")
      id_docalcwog <- "TIMEOUT"
      id_CIDP <- idcase
    }
    if (file.exists((filename))) {
      load(filename)
    } else {
      stop(paste0("ERROR - file not found: ", filename))
    }
    #plot(as(t(saveobj$amat_dag), "graphNEL"))
    #plotAG(saveobj$amat_pag)

    # if (!id_docalcwog) {
    #   L <- saveobj$L
    #   ret_docalcwog <- rundocalcwograph(L, timeout = 1800) #30 minutes of timeout
    #   id_docalcwog <- ret_docalcwog$id
    #   formular <- ret_docalcwog$formular
    #   saveobj$formular <- formular
    #   save(saveobj, file=paste0("./random/ins_", n, "_CIDP_", id_CIDP, "_docalcwog_", id_docalcwog, "_", i, ".RData"))
    # }

    runPAG = FALSE
    if (saveobj$times_CIDP$cipd > 0) {
      print(filename)
      amat <- saveobj$amat_dag
      x <- saveobj$x
      y <- saveobj$y
      z <- saveobj$z
      if (!runPAG) {
        amat_pag <- saveobj$amat_pag
        times_CIDP <- list()
        times_CIDP$pag <- saveobj$times_CIDP$pag
        retCIDP <- runCIDPgivenPAG(amat_pag, x, y, z)
        id_CIDP <- retCIDP$id
        saveobj$times_CIDP$cidp <- retCIDP$time_cidp
        print(saveobj$times_CIDP$cidp)
        #saveobj = list(L=L, ig=ig, amat_dag=amat, amat_pag=amat_pag, x=x, y=y, z=z, times_CIDP=times_CIDP, formular=formular)
        save(saveobj, file=paste0("./random/ins_", n, "_CIDP_", id_CIDP, "_docalcwog_", id_docalcwog, "_", i, ".RData"))
      } else {
        retCIDP <- runPAGandCIDP(amat, x, y, z)
        print(retCIDP$times$cipd)
        print(retCIDP$ret$Qexpr)
        stop("ERROR!!")
      }
    }
  }
}





####################
# Generating plots #
####################

for (idcaseCIDP in c(TRUE, FALSE)) {
  if (idcaseCIDP == TRUE) {
    dirname <- "./random_id_new/"
  } else {
    dirname <- "./random_nonid_new/"
  }

  getdocalctimes <- TRUE
  times_df <- c()
  all_times_docalcwog <- list()
  all_times_cidp <- list()
  for (n in 5:12) {
    for (i in 1:30) {
      idcaseDoCalc <- idcaseCIDP
      #idcaseDoCalc <- FALSE
      filename <- paste0(dirname, "ins_", n, "_CIDP_", idcaseCIDP, "_docalcwog_", idcaseDoCalc, "_", i, ".RData")

      if (!file.exists((filename))) {
        idcaseDoCalc <- "TIMEOUT"
        filename <- paste0(dirname, "ins_", n, "_CIDP_", idcaseCIDP, "_docalcwog_", idcaseDoCalc, "_", i, ".RData")
      }

      if (!file.exists((filename))) {
        idcaseDoCalc <- "PENDING"
        filename <- paste0(dirname, "ins_", n, "_CIDP_", idcaseCIDP, "_docalcwog_", idcaseDoCalc, "_", i, ".RData")
      }

      if (!file.exists((filename))) {
        message(paste0("ERROR - file not found: ", filename))
        next
      }

      load(filename)


      if (getdocalctimes) {
        if (idcaseDoCalc == "PENDING") {
          curtime_pag <- 5000
          curtime_others <- 5000
        } else if (idcaseDoCalc == "TIMEOUT") {
          if (!(length(saveobj$formular$times) == 1 && saveobj$formular$time == "timeout")) {
            cat(paste0("Check TIMEOUT file n=", n, "; i=", i))
          }
          curtime_pag <- 5000 #Inf
          curtime_others <- 5000 #Inf
        } else {
          if (!is.null(saveobj$formular)) {
            curtime_pag <- saveobj$formular$times$fci
            curtime_others <- saveobj$formular$times$total - curtime_pag
          } else {
            cat(paste0("ERROR!!! file n=", n, "; i=", i))
            stop()
          }
        }

        all_times_docalcwog[["pag"]] <- c(all_times_docalcwog[["pag"]], as.numeric(curtime_pag))
        all_times_docalcwog[["others"]] <- c(all_times_docalcwog[["others"]], as.numeric(curtime_others))
      }

      all_times_cidp[["pag"]] <- c(all_times_cidp[["pag"]], as.numeric(saveobj$times_CIDP$pag))
      all_times_cidp[["cidp"]] <- c(all_times_cidp[["cidp"]], as.numeric(saveobj$times_CIDP$cidp))

      times_df <- rbind(times_df, c(n=n, rep=i))
    }
  }

  times_df_cidp_fci <- data.frame(times_df, alg="CIDP", op="FCI", sec=all_times_cidp[["pag"]])
  times_df_cidp_id <- data.frame(times_df, alg="CIDP", op="ID", sec=all_times_cidp[["cidp"]])

  if (getdocalctimes) {
    times_df_hytt_fci <- data.frame(times_df, alg="Hytt2015", op="FCI", sec=all_times_docalcwog[["pag"]])
    times_df_hytt_id <- data.frame(times_df, alg="Hytt2015", op="ID", sec=all_times_docalcwog[["others"]])
    times_df <- rbind(times_df_hytt_fci, times_df_cidp_fci, times_df_hytt_id, times_df_cidp_id)
    times_df$alg <- factor(times_df$alg, levels=c("CIDP", "Hytt2015"))
    times_df$op <- factor(times_df$op, levels=c("FCI", "ID"))
  } else {
    times_df <- rbind(times_df_cidp_fci, times_df_cipd_id)
    times_df$op <- factor(times_df$op, levels=c("FCI", "ID"))
  }
  times_df$n <- factor(times_df$n, levels=5:12)

  if (idcaseCIDP == TRUE) {
    times_df_id <- cbind(times_df, "ID"="ID")
  } else if (idcaseCIDP == FALSE) {
    times_df_nonid <-  cbind(times_df, "ID"="Non-ID")
  }
}

str(times_df_id)
str(times_df_nonid)
times_df <- rbind(times_df_id, times_df_nonid)
times_df$ID <- as.factor(times_df$ID)
str(times_df)

#save(times_df, file="times_df_5to10.RData")
save(times_df, file="times_df_5to12.RData")

load(file="times_df_5to12.RData")
times_df2 <- times_df
times_df <- subset(times_df, op == "ID")

library(scales)
library(ggplot2)
library(dplyr)

aplot <-
  ggplot(aes(y = sec, fill = alg), data = times_df) +
  geom_boxplot(width=3) +
  #ggplot(times_df, aes(y=sec, fill=alg)) +
  #geom_boxplot(width=5) +
  #scale_y_continuous(trans = log_trans(), breaks = base_breaks(),
  #                   labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  ggtitle("Running time by number of variables and identifiability of P(y|do(x),z)") + # bquote( ~ P(y ~ "|" ~ do(x),z))) +
  #facet_wrap(~ n, nrow=1) +
  scale_x_discrete(name = "Number of observed variables", breaks = NULL) +
  theme(legend.position = "bottom", panel.spacing = unit(0.1, "lines"),
        axis.text.y = element_text(size = 12), #, face = "bold"),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16, face = "bold"),
        strip.text = element_text(size = 12),
        legend.title=element_blank(), legend.text=element_text(size=12)) +
  ylab("Time in seconds for 30 randomly generated models") +
  facet_grid(. ~ ID) +
  facet_grid(ID ~ n,  switch = "x") +
  scale_fill_manual(values=c("salmon", "steelblue1"),
                  name="Algorithm",
                  breaks=c("CIDP", "Hytt2015"),
                  labels=c("CIDP", "Hyttinen et al., 2015"))

pdf("CIDPvsHytt5-12.pdf", width = 8.2, height = 5.7)
#pdf("CIDPvsHytt5-10.pdf", width = 7.2, height = 5.5)
print(aplot)
dev.off()

