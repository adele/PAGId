#' @title Algorithm for identifying marginal causal effects from PAGs
#'
#' @description Decides whether the interventional distribution P(y|do(x)) is
#' identifiable or not given a PAG and the observational distribution P(V).
#'
#' @param amat the pcalg adjacency matrix for the PAG.
#' @param x the names for the variables in the treatment set X
#' @param y the names for the variables in the outcome set Y
#' @param verbose a Boolean indicating whether log messages should be printed out.
#'
#' @return a named list where
#'   'id'    - a boolean indicating whether P(y|do(x)) is identifiable or not;
#'   'Qexpr' - list with every required operation in the identification formula as a string;
#'   'Qop'   - list with every required operation in the identification formula;
#'   'query' - a string representing P(y|do(x)). It can be used to access the high-level
#'    identification expression for the query (e.g.,. Qexpr[[query]] and Qop[[query]])
#'
#' @seealso \code{\link{CIDP}}
#'
#' @export IDP
IDP <- function(amat, x, y, verbose=FALSE) {

  identify <- function(c, t) {
    if (verbose)
      print(paste0("#### Computing Q[", paste0(c, collapse=","), "], from Q[", paste0(t, collapse=","), "]"))

    vnames <- colnames(amat)
    intervset <- setdiff(vnames, c)
    intervset_str <- "obs"
    Pinterv <- "P"
    if( length(intervset) > 0) {
      intervset_str <- paste0(intervset, collapse=",")
      Pinterv <- paste0("P_{", intervset_str, "}")
    }

    if (length(c) == 0) {
      if (verbose)
        print("Returning 1")
      QExprList[[intervset_str]] <<- 1
      QOpList[[intervset_str]] <<- 1
      return()
    }

    if (length(setdiff(t, c)) == 0) {
      if (verbose)
        print("Returning Q[C] = Q[T]")
      if (is.null(QExprList[[intervset_str]])) {
        QExprList[[intervset_str]] <<- paste0(Pinterv, "(", paste0(t, collapse=","), ")")
        QOpList[[intervset_str]] <<- list(type="none", param=list(interv=intervset, var=t))
      }
      return()
    }

    amat_t <- inducedPag(amat, t, F)$amat
    amat_c <- inducedPag(amat, c, F)$amat
    tminusc <- setdiff(t, c)
    buckets_t <- getBucketList(amat_t)

    for (b_t in buckets_t) {
      if (all(b_t %in% tminusc)) { # B is contained in T\C
        pcb <- getPCComponentA(amat, amat_t, b_t) # pc-component of the bucket B in P_T
        possDeB <- getPossDescendants(amat_t, b_t) # possible descendants of the bucket B in P_T
        interB <- intersect(pcb, possDeB)

        if (all(interB %in% b_t)) { # interB contained in B
          if (verbose)
            print(paste0("Cond. in line 6 is satisfied for bucket {", paste0(b_t, collapse=","), "}"))

          tminusb <- setdiff(t, b_t)
          intervset <- setdiff(vnames, tminusb)
          intervset_str <- "obs"
          Pinterv <- "P"
          if( length(intervset) > 0) {
            intervset_str <- paste0(intervset, collapse=",")
            Pinterv <- paste0("P_{", intervset_str, "}")
          }
          if (verbose)
            print(paste0("Q[T\\B] = ", Pinterv))

          vminust <- setdiff(vnames, t)
          P_vminust <- "P" # Q[T}
          if (length(vminust) > 0) {
            P_vminust <- paste0("P_{", paste0(vminust, collapse=","), "}")
          }
          if (verbose)
            print(paste0("Computing Q[T\\B]=Q[", paste0(tminusb, collapse=","), "]=", Pinterv,
                         " from Q=Q[T]=Q[", paste0(t, collapse=","), "]=", P_vminust,
                         " via Prop. 6"))

          # Q[T\\B] = P_{vminust}(t) / P_{vminust}(b_t | t\possDeB )
          QOpList[[intervset_str]] <<-
            list(type="prop6", param=list(interv=vminust, num.var=t, den.var=b_t, den.cond=setdiff(t, possDeB)))

          Qtminusb_str <-
            paste0("\\frac{", P_vminust, "(", paste0(t, collapse=","), ")","}{",
                   P_vminust, "(", paste0(b_t, collapse=","), " | ", paste0(setdiff(t, possDeB), collapse = ","), ")}")

          QExprList[[intervset_str]] <<- Qtminusb_str

          identify(c, tminusb)

          return()
        }
      }
    }

    for (b_t in buckets_t) {
      if (all(b_t %in% c)) { # B is contained in C
        region_bt <- getRegion(amat, amat_c, b_t)
        if (length(setdiff(c, region_bt)) != 0) { # Region of B != C
          if (verbose)
            print(paste0("Cond. in line 9 is satisfied for bucket {", paste0(b_t, collapse=","), "}"))

          region_cminusrb <- getRegion(amat, amat_c, setdiff(c, region_bt))
          region_inter <- intersect(region_bt, region_cminusrb)

          if (verbose) {
            print(paste0("Computing Q[C]=Q[", paste0(c, collapse=","), "]=", Pinterv,
                         " from via Prop. 7"))
            print(paste0("R_b: ", paste0(region_bt, collapse=",")))
            print(paste0("R_cminusrb: ", paste0(region_cminusrb, collapse=",")))
            print(paste0("R_inter: ", paste0(region_inter, collapse=",")))
          }

          q_rb <- identify(region_bt, t)
          q_cminusrb <- identify(region_cminusrb, t)
          q_inter <- identify(region_inter, t)

          rb_interv <- setdiff(vnames, region_bt)
          rcmb_interv <- setdiff(vnames, region_cminusrb)
          rinter_interv <- setdiff(vnames, region_inter)

          # Q[C] = (Q[Rb] * Q[R_{C\Rb}]) / Q[Rinter]

          QOpList[[intervset_str]] <<- list(type="prop7", param=list(num.prod1=rb_interv, num.prod2=rcmb_interv, den=rinter_interv))

          QRb_str = Qcminusrb_str = Qc_str = "P"
          if (length(rb_interv) > 0) {
            QRb_str <- paste0("P_{", paste0(rb_interv, collapse = ","), "}")
          }
          if (length(rcmb_interv) > 0) {
            Qcminusrb_str <- paste0("P_{", paste0(rcmb_interv, collapse = ","), "}")
          }
          if (length(rinter_interv) > 0) {
            Qinter_str <- paste0("P_{", paste0(rinter_interv, collapse = ","), "}")
          }
          Qc_str <- paste0("\\frac{", QRb_str, " . ",  Qcminusrb_str  ,"}{", Qinter_str, "}")

          QExprList[[intervset_str]] <<- Qc_str

          return()
        }
      }
    }

    stop(paste0("Q[", paste0(c, collapse=","), "] is not identifiable from Q[", paste0(t, collapse=","), "]"))
  }

  QExprList <- list()
  QOpList <- list()

  v <- colnames(amat)
  vminusx <- setdiff(v, x)
  amat_vminusx <- inducedPag(amat, vminusx)$amat

  d <- getPossAncestors(amat_vminusx, y)
  dminusy <- setdiff(d, y)

  id <- tryCatch({
    identify(d, v) # Q[d] = P_{v\d}
    TRUE
    },
    error = function(err) {
      print(paste("FAIL: ", err))
      FALSE
    }
  )

  intervset_str <- "obs"
  Pinterv <- "P"
  if( length(x) > 0) {
    intervset_str <- paste0(x, collapse=",")
    Pinterv <- paste0("P_{", intervset_str, "}")
  }
  query <- paste0(Pinterv, "(", paste0(y, collapse=","), ")")


  qd_intervset <- setdiff(v, d) # Q[d] = P_{v\d}
  qd_intervset_str <- "obs"
  if (length(qd_intervset) > 0) {
    qd_intervset_str <- paste0(qd_intervset, collapse=",")
  }

  if (id) {
    sumset_exp <- ""
    if (length(dminusy) > 0) {
      sumset_exp <- paste0("\\sum_{", paste0(dminusy, collapse=","), "}")
      QExprList[[query]] <- paste0(sumset_exp, QExprList[[qd_intervset_str]])
      QOpList[[query]] <- list(type="sumset", param=list(sumset=dminusy, interv=qd_intervset))
    } else {
      QExprList[[query]] <- QExprList[[qd_intervset_str]]
      QExprList <- QExprList[-which(names(QExprList) == qd_intervset_str)]
      QOpList[[query]] <- QOpList[[qd_intervset_str]]
      QOpList <- QOpList[-which(names(QOpList) == qd_intervset_str)]
    }
    return(list(id=TRUE, query=query, Qop=QOpList, Qexpr=QExprList))
  } else {
    return(list(id=FALSE, query=query))
  }
}

