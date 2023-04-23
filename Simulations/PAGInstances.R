# For PAGs, the edgemark-code refers to the column index
# 0: No edge
# 1: Circle
# 2: Arrowhead
# 3: Tail

getDAG2aHyttinen <- function(ret_dagg = TRUE, simpleDAG = FALSE) {
  allvars <- c("x", "y", "z", "h", "w", "uwx", "uxh")
  p <- length(allvars)
  amat <- matrix(0, p, p)
  colnames(amat) <- rownames(amat) <- allvars
  amat["z","y"] <- 0; amat["y","z"] <- 1; # z -> y
  amat["x","z"] <- 0; amat["z","x"] <- 1; # x -> z
  amat["h","y"] <- 0; amat["y","h"] <- 1; # h -> y
  amat["uwx","x"] <- 0; amat["x","uwx"] <- 1; # uwx -> x
  amat["uwx","w"] <- 0; amat["w","uwx"] <- 1; # uwx -> w
  amat["uxh","h"] <- 0; amat["h","uxh"] <- 1; # uxh -> h
  amat["uxh","x"] <- 0; amat["x","uxh"] <- 1; # uxh -> x
  lat <- c("uwx", "uxh")

  #plot(as(t(amat), "graphNEL"))

  x <- c("x")
  y <- c("y")
  z <- c()

  dagg=NULL
  if (ret_dagg) {
    dagg <- pcalg::pcalg2dagitty(amat, colnames(amat), type="dag")
    #plot(dagg)

    dagitty::latents(dagg) <- lat
    dagitty::exposures(dagg) <- x
    dagitty::outcomes(dagg) <- y
  }
  #plot(getPAG(dagg)) # same in Fig. 5c in Hyttinen et al., 2015

  return(list(amat=amat, x=x, y=y, z=z, lat=lat, dagg=dagg))
}


getDAG1 <- function(ret_dagg = TRUE, simpleDAG = FALSE) {
  if (simpleDAG) {
    allvars <- c("x1", "x2", "y", "a", "w", "z")
    p <- length(allvars)
    amat <- matrix(0, p, p)
    colnames(amat) <- rownames(amat) <- allvars
    amat["w","z"] <- 0; amat["z","w"] <- 1; # w -> z
    amat["z","y"] <- 0; amat["y","z"] <- 1; # z -> y
    amat["a","y"] <- 0; amat["y","a"] <- 1; # a -> y
    amat["a","w"] <- 0; amat["w","a"] <- 1; # a -> w
    amat["a","x1"] <- 0; amat["x1","a"] <- 1; # a -> x1
    amat["x2","x1"] <- 0; amat["x1","x2"] <- 1; # x2 -> x1
    amat["x2","w"] <- 0; amat["w","x2"] <- 1; # x2 -> w
    amat["x1","w"] <- 0; amat["w","x1"] <- 1; # x1 -> w
    lat <- c()
  } else {
    allvars <- c("x1", "x2", "y", "a", "w", "z", "uaw", "uax1", "uwx2", "ux1x2")
    p <- length(allvars)
    amat <- matrix(0, p, p)
    colnames(amat) <- rownames(amat) <- allvars
    amat["w","z"] <- 0; amat["z","w"] <- 1; # w -> z
    amat["z","y"] <- 0; amat["y","z"] <- 1; # z -> y
    amat["a","y"] <- 0; amat["y","a"] <- 1; # a -> y
    amat["a","w"] <- 0; amat["w","a"] <- 1; # a -> w
    amat["a","x1"] <- 0; amat["x1","a"] <- 1; # a -> x1
    amat["x1","w"] <- 0; amat["w","x1"] <- 1; # x1 -> w
    amat["x2","x1"] <- 0; amat["x1","x2"] <- 1; # x2 -> x1
    amat["ux1x2","x1"] <- 0; amat["x1","ux1x2"] <- 1; # ux1x2 -> x1
    amat["ux1x2","x2"] <- 0; amat["x2","ux1x2"] <- 1; # ux1x2 -> x2
    amat["uwx2","w"] <- 0; amat["w","uwx2"] <- 1; # uwx2 -> w
    amat["uwx2","x2"] <- 0; amat["x2","uwx2"] <- 1; # uwx2 -> x2
    amat["uax1","a"] <- 0; amat["a","uax1"] <- 1; # uax1 -> a
    amat["uax1","x1"] <- 0; amat["x1","uax1"] <- 1; # uax1 -> x1
    amat["uaw","a"] <- 0; amat["a","uaw"] <- 1; # uaw -> a
    amat["uaw","w"] <- 0; amat["w","uaw"] <- 1; # uaw -> w
    lat <- c("uaw", "uax1", "uwx2", "ux1x2")
  }
  #plot(as(t(amat), "graphNEL"))

  x <- c("x1", "x2", "w")
  y <- c("y")
  z <- c("z")

  dagg=NULL
  if (ret_dagg) {
   dagg <- pcalg::pcalg2dagitty(amat, colnames(amat), type="dag")
   #plot(dagg)

   dagitty::latents(dagg) <- lat
   dagitty::exposures(dagg) <- x
   dagitty::outcomes(dagg) <- y
  }
  #plot(getPAG(dagg))

  return(list(amat=amat, x=x, y=y, z=z, lat=lat, dagg=dagg))
}


getPAG1 <- function(ret_dagg = TRUE) {
  allvars <- c("x1", "x2", "y", "a", "w", "z")
  p <- length(allvars)
  amat <- matrix(0, p, p)
  colnames(amat) <- rownames(amat) <- allvars
  amat["w","z"] <- 2; amat["z","w"] <- 3; # w -> z
  amat["z","y"] <- 2; amat["y","z"] <- 3; # z -> y
  amat["a","y"] <- 2; amat["y","a"] <- 3; # a -> y
  amat["a","w"] <- 2; amat["w","a"] <- 1; # a @-> w
  amat["a","x1"] <- 2; amat["x1","a"] <- 1; # a @-> x1
  amat["x2","x1"] <- 2; amat["x1","x2"] <- 1; # x2 @-> x1
  amat["x2","w"] <- 2; amat["w","x2"] <- 1; # x2 @-> w
  amat["x1","w"] <- 1; amat["w","x1"] <- 1; # x1 @-> w
  #plotAG(amat)

  x <- c("w", "x1", "x2")
  y <- c("y")
  z <- c("z")
  return(list(amat=amat, x=x, y=y, z=z))
}


getDAG2n <- function(ret_dagg = TRUE, simpleDAG = FALSE) {
  if (simpleDAG) {
    allvars <- c("a", "b", "c", "d", "e", "f", "y")
    p <- length(allvars)
    amat <- matrix(0, p, p)
    colnames(amat) <- rownames(amat) <- allvars
    amat["e","y"] <- 0; amat["y","e"] <- 1; # e -> y
    amat["b","a"] <- 0; amat["a","b"] <- 1; # b -> a
    amat["b","y"] <- 0; amat["y","b"] <- 1; # b -> y
    amat["b","f"] <- 0; amat["f","b"] <- 1; # b -> f
    amat["b","e"] <- 0; amat["e","b"] <- 1; # b -> e
    amat["f","e"] <- 0; amat["e","f"] <- 1; # e -> f
    amat["f","y"] <- 0; amat["y","f"] <- 1; # f -> y
    amat["d","e"] <- 0; amat["e","d"] <- 1; # d -> e
    amat["d","f"] <- 0; amat["f","d"] <- 1; # d -> f
    amat["c","e"] <- 0; amat["e","c"] <- 1; # c -> e
    lat <- c()
  } else {
    allvars <- c("a", "b", "c", "d", "e", "f", "y", "uab")
    p <- length(allvars)
    amat <- matrix(0, p, p)
    colnames(amat) <- rownames(amat) <- allvars
    amat["e","y"] <- 0; amat["y","e"] <- 1; # e -> y
    amat["b","y"] <- 0; amat["y","b"] <- 1; # b -> y
    amat["b","f"] <- 0; amat["f","b"] <- 1; # b -> f
    amat["b","e"] <- 0; amat["e","b"] <- 1; # b -> e
    amat["f","e"] <- 0; amat["e","f"] <- 1; # e -> f
    amat["f","y"] <- 0; amat["y","f"] <- 1; # f -> y
    amat["d","e"] <- 0; amat["e","d"] <- 1; # d -> e
    amat["d","f"] <- 0; amat["f","d"] <- 1; # d -> f
    amat["c","e"] <- 0; amat["e","c"] <- 1; # c -> e
    amat["uab","a"] <- 0; amat["a","uab"] <- 1; # uab -> a
    amat["uab","b"] <- 0; amat["b","uab"] <- 1; # uab -> b
    lat <- c("uab")
  }
  #plot(as(t(amat), "graphNEL"))

  #P_{a,f}(y|b,e)
  x <- c("a", "f")
  y <- c("y")
  z <- c("b", "e")

  dagg=NULL
  if (ret_dagg) {
    dagg <- pcalg::pcalg2dagitty(amat, colnames(amat), type="dag")
    #plot(dagg)

    dagitty::latents(dagg) <- lat
    dagitty::exposures(dagg) <- x
    dagitty::outcomes(dagg) <- y
  }
  #plot(getPAG(dagg))

  return(list(amat=amat, x=x, y=y, z=z, lat=lat, dagg=dagg))
}


getDAG2s <- function(ret_dagg = TRUE, simpleDAG = FALSE) {
  if (simpleDAG) {
    allvars <- c("a", "b", "c", "e", "y")
    p <- length(allvars)
    amat <- matrix(0, p, p)
    colnames(amat) <- rownames(amat) <- allvars
    amat["b","y"] <- 0; amat["y","b"] <- 1; # b -> y
    amat["e","y"] <- 0; amat["y","e"] <- 1; # e -> y
    amat["c","y"] <- 0; amat["y","c"] <- 1; # c -> y
    amat["a","e"] <- 0; amat["e","a"] <- 1; # a -> e
    amat["a","c"] <- 0; amat["c","a"] <- 1; # a -> c
    amat["b","e"] <- 0; amat["e","b"] <- 1; # b -> e
    amat["b","c"] <- 0; amat["c","b"] <- 1; # b -> c
    amat["c","e"] <- 0; amat["e","c"] <- 1; # c -> e
    lat <- c()
  } else {
    allvars <- c("a", "b", "c", "e", "y", "uce", "ucb", "uae")
    p <- length(allvars)
    amat <- matrix(0, p, p)
    colnames(amat) <- rownames(amat) <- allvars
    amat["a","e"] <- 0; amat["e","a"] <- 1; # a -> e
    amat["a","c"] <- 0; amat["c","a"] <- 1; # a -> c
    amat["b","e"] <- 0; amat["e","b"] <- 1; # b -> e
    amat["b","c"] <- 0; amat["c","b"] <- 1; # b -> c
    amat["b","y"] <- 0; amat["y","b"] <- 1; # b -> y
    amat["c","e"] <- 0; amat["e","c"] <- 1; # c -> e
    amat["c","y"] <- 0; amat["y","c"] <- 1; # c -> y
    amat["e","y"] <- 0; amat["y","e"] <- 1; # e -> y
    amat["uce","c"] <- 0; amat["c","uce"] <- 1; # uce -> c
    amat["uce","e"] <- 0; amat["e","uce"] <- 1; # uce -> e
    amat["ucb","c"] <- 0; amat["c","ucb"] <- 1; # ucb -> c
    amat["ucb","b"] <- 0; amat["b","ucb"] <- 1; # ucb -> b
    amat["uae","a"] <- 0; amat["a","uae"] <- 1; # uae -> a
    amat["uae","e"] <- 0; amat["e","uae"] <- 1; # uae -> e
    lat <- c("uce", "ucb", "uae")
  }
  #plot(as(t(amat), "graphNEL"))

  x <- c("c")
  y <- c("y")
  z <- c("b", "e")

  dagg=NULL
  if (ret_dagg) {
    dagg <- pcalg::pcalg2dagitty(amat, colnames(amat), type="dag")
    #plot(dagg)

    dagitty::latents(dagg) <- lat
    dagitty::exposures(dagg) <- x
    dagitty::outcomes(dagg) <- y
  }
  #plot(getPAG(dagg))

  return(list(amat=amat, x=x, y=y, z=z, lat=lat, dagg=dagg))
}


getPAG2s <- function(ret_dagg = TRUE) {
  allvars <- c("a", "b", "c", "e", "y")
  p <- length(allvars)
  amat <- matrix(0, p, p)
  colnames(amat) <- rownames(amat) <- allvars
  amat["b","y"] <- 2; amat["y","b"] <- 3; # z -> y
  amat["e","y"] <- 2; amat["y","e"] <- 3; # w -> y
  amat["c","y"] <- 2; amat["y","c"] <- 3; # w -> y
  amat["a","e"] <- 2; amat["e","a"] <- 1; # a @-> e
  amat["a","c"] <- 2; amat["c","a"] <- 1; # a @-> c
  amat["b","e"] <- 2; amat["e","b"] <- 1; # b @-> e
  amat["b","c"] <- 2; amat["c","b"] <- 1; # b @-> c
  amat["c","e"] <- 1; amat["e","c"] <- 1; # c @-@ e
  #plotAG(amat)

  x <- c("c")
  y <- c("y")
  z <- c("b", "e")
  return(list(amat=amat, x=x, y=y, z=z))
}


getPAG3 <- function(ret_dagg = TRUE) {
  allvars <- c("x", "w", "y", "z", "t")
  p <- length(allvars)
  amat <- matrix(0, p, p)
  colnames(amat) <- rownames(amat) <- allvars
  amat["z","y"] <- 2; amat["y","z"] <- 3; # z -> y
  amat["w","y"] <- 2; amat["y","w"] <- 3; # w -> y
  amat["w","x"] <- 2; amat["x","w"] <- 1; # w @-> x
  amat["w","z"] <- 2; amat["z","w"] <- 1; # w @-> z
  amat["t","x"] <- 2; amat["x","t"] <- 1; # t @-> x
  amat["t","z"] <- 2; amat["z","t"] <- 1; # t @-> z
  amat["x","z"] <- 1; amat["z","x"] <- 1; # x @-@ z
  #plotAG(amat)

  x <- c("x")
  y <- c("y")
  z <- c("z")
  return(list(amat=amat, x=x, y=y, z=z))
}

# DAG in Fig. 2b), represented by PAG in Fig. 2a) in R50 - https://causalai.net/r50.pdf
getDAG2bR42 <- function(ret_dagg = TRUE) {
  allvars <- c("x1", "x2", "y1", "y2", "y3", "y4", "y5", "ux1x2", "uy2y3", "uy4y5", "ux1y3", "ux2y1")
  p <- length(allvars)
  amat <- matrix(0, p, p)
  colnames(amat) <- rownames(amat) <- allvars
  amat["x2", "y2"] <- 0; amat["y2", "x2"] <- 1; # x2 -> y2
  amat["x1", "y1"] <- 0; amat["y1", "x1"] <- 1; # x1 -> y1
  amat["y4", "y3"] <- 0; amat["y3", "y4"] <- 1; # y4 -> y3
  amat["y5", "y1"] <- 0; amat["y1", "y5"] <- 1; # y5 -> y1
  amat["y5", "y4"] <- 0; amat["y4", "y5"] <- 1; # y5 -> y4
  amat["ux1x2", "x1"] <- 0; amat["x1", "ux1x2"] <- 1; # ux1x2 -> x1
  amat["ux1x2", "x2"] <- 0; amat["x2", "ux1x2"] <- 1; # ux1x2 -> x2
  amat["uy2y3", "y2"] <- 0; amat["y2", "uy2y3"] <- 1; # uy2y3 -> y2
  amat["uy2y3", "y3"] <- 0; amat["y3", "uy2y3"] <- 1; # uy2y3 -> y3
  amat["uy4y5", "y4"] <- 0; amat["y4", "uy4y5"] <- 1; # uy4y5 -> y4
  amat["uy4y5", "y5"] <- 0; amat["y5", "uy4y5"] <- 1; # uy4y5 -> y5
  amat["ux1y3", "x1"] <- 0; amat["x1", "ux1y3"] <- 1; # ux1y3 -> x1
  amat["ux1y3", "y3"] <- 0; amat["y3", "ux1y3"] <- 1; # ux1y3 -> y3
  amat["ux2y1", "x2"] <- 0; amat["x2", "ux2y1"] <- 1; # ux2y1 -> x2
  amat["ux2y1", "y1"] <- 0; amat["y1", "ux2y1"] <- 1; # ux2y1 -> y1
  # plot(as(t(amat), "graphNEL"))

  x <- c("x1", "x2")
  y <- c("y1", "y2", "y3", "y4", "y5")
  z <- NULL
  lat <- c("ux1x2", "uy2y3", "uy4y5", "ux1y3", "ux2y1")

  dagg=NULL
  if (ret_dagg) {
    dagg <- pcalg::pcalg2dagitty(amat, colnames(amat),type="dag")
    #plot(dagg)

    dagitty::latents(dagg) <- lat
    dagitty::exposures(dagg) <- x
    dagitty::outcomes(dagg) <- y
  }
  #plot(getPAG(dagg))

  return(list(amat=amat, x=x, y=y, z=z, lat=lat, dagg=dagg))
}


# PAG in Fig. 2a) in R42 - https://causalai.net/r42.pdf
getPAG2aR42 <- function() {
  amat <- matrix(0, 7, 7)
  colnames(amat) <- c("x1", "x2", "y1", "y2", "y3", "y4", "y5")
  rownames(amat) <- c("x1", "x2", "y1", "y2", "y3", "y4", "y5")
  amat["x1","x2"] <- 2; amat["x2","x1"] <- 2; #x1 <-> x2
  amat["x1","y3"] <- 2; amat["y3","x1"] <- 2; #x1 <-> y3
  amat["x1","y1"] <- 2; amat["y1","x1"] <- 3; #x1 -> y1
  amat["x2","y2"] <- 2; amat["y2","x2"] <- 3; #x2 -> y2
  amat["x2","y1"] <- 2; amat["y1","x2"] <- 2; #x2 <-> y1
  amat["y2","y3"] <- 2; amat["y3","y2"] <- 2; #y2 <-> y3
  amat["y3","y4"] <- 1; amat["y4","y3"] <- 2; #y3 <-@ y4
  amat["y4","y5"] <- 1; amat["y5","y4"] <- 1; #y4 @-@ y5
  amat["y1","y5"] <- 3; amat["y5","y1"] <- 2; #y1 <- y5
  #plotAG(amat)

  x <- c("x1", "x2")
  y <- c("y1", "y2", "y3", "y4", "y5")
  z <- c()
  return(list(amat=amat, x=x, y=y, z=z))
}


# PAG in Fig. 3a) in R42 - https://causalai.net/r42.pdf
# Px(y) is not identifiable
getPAG3aR42 <- function() {
  amat <- matrix(0, 2, 2)
  colnames(amat) <- c("x", "y")
  rownames(amat) <- c("x", "y")
  amat["x","y"] <- 1; amat["y","x"] <- 1; # x o-o y
  #plotAG(amat)

  x <- c("x")
  y <- c("y")
  z <- c("")
  return(list(amat=amat, x=x, y=y, z=z))
}

# PAG in Fig. 3b) in R42 - https://causalai.net/r42.pdf
# Px(y) is not identifiable
getPAG3bR42 <- function() {
  amat <- matrix(0, 4, 4)
  colnames(amat) <- c("x", "v1", "y1", "y2")
  rownames(amat) <- c("x", "v1", "y1", "y2")
  amat["x","y2"] <- 2; amat["y2","x"] <- 3; # x -> y2
  amat["x","y1"] <- 2; amat["y1","x"] <- 2; # x <-> y1
  amat["y2","y1"] <- 2; amat["y1","y2"] <- 2; # y2 <-> y1
  amat["x","v1"] <- 1; amat["v1","x"] <- 2; # v1 o-> x
  #plotAG(amat)

  x <- c("x")
  y <- c("y1", "y2")
  z <- c("")
  return(list(amat=amat, x=x, y=y, z=z))
}


# PAG in Fig. 3c) in R42 - https://causalai.net/r42.pdf
# Px(y) is not identifiable
getPAG3cR42 <- function() {
  amat <- matrix(0, 6, 6)
  colnames(amat) <- c("x1", "x2", "y1", "y2", "y3", "y4")
  rownames(amat) <- c("x1", "x2", "y1", "y2", "y3", "y4")
  amat["x1","x2"] <- 2; amat["x2","x1"] <- 2; #x1 <-> x2
  amat["x1","y3"] <- 2; amat["y3","x1"] <- 2; #x1 <-> y3
  amat["x1","y1"] <- 2; amat["y1","x1"] <- 3; #x1 -> y1
  amat["x2","y2"] <- 2; amat["y2","x2"] <- 3; #x2 -> y2
  amat["x2","y1"] <- 2; amat["y1","x2"] <- 2; #x2 <-> y1
  amat["y2","y3"] <- 2; amat["y3","y2"] <- 2; #y2 <-> y3
  amat["y3","y4"] <- 2; amat["y4","y3"] <- 2; #y3 <-> y4
  amat["y1","y4"] <- 2; amat["y4","y1"] <- 2; #y1 <-> y4
  #plotAG(amat)

  x <- c("x1", "x2")
  y <- c("y1", "y2", "y3", "y4")
  z <- c()
  return(list(amat=amat, x=x, y=y, z=z))
}

# A DAG represented by the PAG in Fig. 2a) in R50 (Jaber et al., 2019)
getDAG2aR50 <- function(ret_dagg = TRUE) {
  allvars <- c("x", "y", "z", "uzy")
  p <- length(allvars)
  amat <- matrix(0,p,p) # notation: column causes row
  rownames(amat) <- colnames(amat) <- allvars
  amat["z", "x"] <- 0; amat["x", "z"] <- 1; # z -> x
  amat["z", "y"] <- 0; amat["y", "z"] <- 1; # z -> y
  amat["uzy", "z"] <- 0; amat["z", "uzy"] <- 1; # uzy -> z
  amat["uzy", "y"] <- 0; amat["y", "uzy"] <- 1; # uzy -> y
  # plot(as(t(amat), "graphNEL"))

  x <- c("x")
  y <- c("y")
  z <- c("z")
  lat <- c("uzy")

  dagg=NULL
  if (ret_dagg) {
    dagg <- pcalg::pcalg2dagitty(amat, colnames(amat),type="dag")
    #plot(dagg)

    dagitty::latents(dagg) <- lat
    dagitty::exposures(dagg) <- x
    dagitty::outcomes(dagg) <- y
  }
  #plot(getPAG(dagg))

  return(list(amat=amat, x=x, y=y, z=z, lat=lat, dagg=dagg))
}


# PAG in Fig. 2a) in R50 - https://causalai.net/r50.pdf
getPAG2aR50 <- function() {
  amat <- matrix(0, 3, 3)
  colnames(amat) <- c("x", "z", "y")
  rownames(amat) <- c("x", "z", "y")
  amat["x","z"] <- 1; amat["z","x"] <- 1; # x o-o z
  amat["y","z"] <- 1; amat["z","y"] <- 1; # y o-o z
  #plotAG(amat)

  x <- c("x")
  y <- c("y")
  z <- c("z")
  return(list(amat=amat, x=x, y=y, z=z))
}


# A DAG represented by the PAG in Fig. 2b) in R50 (Jaber et al., 2019)
getDAG2bR50 <- function(ret_dagg = TRUE) {
  allvars <- c("x", "y", "z", "w", "uzw", "uzy")
  p <- length(allvars)
  amat <- matrix(0,p,p) # notation: column causes row
  rownames(amat) <- colnames(amat) <- allvars
  amat["z", "x"] <- 0; amat["x", "z"] <- 1; # z -> x
  amat["z", "y"] <- 0; amat["y", "z"] <- 1; # z -> y
  amat["w", "x"] <- 0; amat["x", "w"] <- 1; # w -> x
  amat["w", "y"] <- 0; amat["y", "w"] <- 1; # w -> y
  amat["uzw", "z"] <- 0; amat["z", "uzw"] <- 1; # uzw -> z
  amat["uzw", "w"] <- 0; amat["w", "uzw"] <- 1; # uzw -> w
  amat["uzy", "z"] <- 0; amat["z", "uzy"] <- 1; # uzy -> z
  amat["uzy", "y"] <- 0; amat["y", "uzy"] <- 1; # uzy -> y
  # plot(as(t(amat), "graphNEL"))

  x <- c("x")
  y <- c("y")
  z <- c("z", "w")
  lat <- c("uzw", "uzy")

  dagg=NULL
  if (ret_dagg) {
    dagg <- pcalg::pcalg2dagitty(amat, colnames(amat),type="dag")
    #plot(dagg)

    dagitty::latents(dagg) <- lat
    dagitty::exposures(dagg) <- x
    dagitty::outcomes(dagg) <- y
  }
  #plot(getPAG(dagg))

  return(list(amat=amat, x=x, y=y, z=z, lat=lat, dagg=dagg))
}

# PAG in Fig. 2b) in R50 - https://causalai.net/r50.pdf
getPAG2bR50 <- function() {
  amat <- matrix(0, 4, 4)
  colnames(amat) <- c("x", "z", "w", "y")
  rownames(amat) <- c("x", "z", "w", "y")
  amat["x","z"] <- 1; amat["z","x"] <- 1; # x o-o z
  amat["w","z"] <- 1; amat["z","w"] <- 1; # w o-o z
  amat["y","z"] <- 1; amat["z","y"] <- 1; # y o-o z
  amat["x","w"] <- 1; amat["w","x"] <- 1; # x o-o w
  amat["y","w"] <- 1; amat["w","y"] <- 1; # y o-o w
  #plotAG(amat)

  x <- c("x")
  y <- c("y")
  z <- c("z", "w")
  return(list(amat=amat, x=x, y=y, z=z))
}

# A DAG represented by the PAG in Fig. 2c) in R50 - https://causalai.net/r50.pdf
getDAG2cR50 <- function(ret_dagg = TRUE) {
  p <- 6
  amat <- matrix(0,p,p) # notation: column causes row
  rownames(amat) <- colnames(amat) <- c("x", "z1", "z2", "y", "u1", "u2")
  amat["u2", "z1"] <- 0; amat["z1", "u2"] <- 1; # u2 -> z1
  amat["u2", "x"] <- 0; amat["x", "u2"] <- 1; # u2 -> x
  amat["x", "z1"] <- 0; amat["z1", "x"] <- 1; # x -> z1
  amat["u1", "z1"] <- 0; amat["z1", "u1"] <- 1; # u1 -> z1
  amat["u1", "z2"] <- 0; amat["z2", "u1"] <- 1; # u1 -> z2
  amat["z2", "z1"] <- 0; amat["z1", "z2"] <- 1; # z2 -> z1
  amat["z2", "y"] <- 0; amat["y", "z2"] <- 1; # z2 -> y
  amat["y", "z1"] <- 1; amat["z1", "y"] <- 0; # z1 -> y
  # plot(as(t(amat), "graphNEL"))

  x <- c("x")
  y <- c("y")
  z <- c("z1", "z2")
  lat <- c("u1", "u2")

  dagg=NULL
  if (ret_dagg) {
    dagg <- pcalg::pcalg2dagitty(amat, colnames(amat),type="dag")
    #plot(dagg)

    dagitty::latents(dagg) <- lat
    dagitty::exposures(dagg) <- x
    dagitty::outcomes(dagg) <- y
  }
  #plot(getPAG(dagg))

  return(list(amat=amat, x=x, y=y, z=z, lat=lat, dagg=dagg))
}


# PAG in Fig. 2c) in R50 - https://causalai.net/r50.pdf
getPAG2cR50 <- function() {
  amat <- matrix(0, 4, 4)
  colnames(amat) <- c("x", "z1", "z2", "y")
  rownames(amat) <- c("x", "z1", "z2", "y")
  amat["x","z1"] <- 2; amat["z1","x"] <- 1; # x o-> z1
  amat["z2","z1"] <- 2; amat["z1","z2"] <- 1; # z2 o-> z1
  amat["z1","y"] <- 2; amat["y","z1"] <- 3; # z1 -> y
  amat["z2","y"] <- 2; amat["y","z2"] <- 3; # z2 -> y
  #plotAG(amat)

  x <- c("x")
  y <- c("y")
  z <- c("z1", "z2")
  return(list(amat=amat, x=x, y=y, z=z))
}


# A DAG represented by the PAG in Fig. 2c) in R50 - https://causalai.net/r50.pdf
getDAG2dR50 <- function(ret_dagg = TRUE) {
  allvars <- c("x", "z1", "z2", "w", "y", "uxz1", "uwz2")
  p <- length(allvars)
  amat <- matrix(0,p,p) # notation: column causes row
  rownames(amat) <- colnames(amat) <- allvars
  amat["uxz1", "z1"] <- 0; amat["z1", "uxz1"] <- 1; # uxz1 -> z1
  amat["uxz1", "x"] <- 0; amat["x", "uxz1"] <- 1; # uxz1 -> x
  amat["x", "z1"] <- 0; amat["z1", "x"] <- 1; # x -> z1
  amat["w", "z1"] <- 0; amat["z1", "w"] <- 1; # w -> z1
  amat["uwz2", "z2"] <- 0; amat["z2", "uwz2"] <- 1; # uwz2 -> z2
  amat["uwz2", "w"] <- 0; amat["w", "uwz2"] <- 1; # uwz2 -> w
  amat["z2", "y"] <- 0; amat["y", "z2"] <- 1; # z2 -> y
  amat["y", "z1"] <- 1; amat["z1", "y"] <- 0; # z1 -> y
  # plot(as(t(amat), "graphNEL"))

  x <- c("x")
  y <- c("y")
  z <- c("z1", "z2")
  lat <- c("uwz2", "uxz1")

  dagg=NULL
  if (ret_dagg) {
    dagg <- pcalg::pcalg2dagitty(amat, colnames(amat),type="dag")
    #plot(dagg)

    dagitty::latents(dagg) <- lat
    dagitty::exposures(dagg) <- x
    dagitty::outcomes(dagg) <- y
  }
  #plot(getPAG(dagg))

  return(list(amat=amat, x=x, y=y, z=z, lat=lat, dagg=dagg))
}

# PAG in Fig. 2d) in R50 - https://causalai.net/r50.pdf
getPAG2dR50 <- function() {
  amat <- matrix(0, 5, 5)
  colnames(amat) <- c("x", "z1", "w", "z2", "y")
  rownames(amat) <- c("x", "z1", "w", "z2", "y")
  amat["x","z1"] <- 2; amat["z1","x"] <- 1; # x o-> z1
  amat["w","z1"] <- 2; amat["z1","w"] <- 1; # w o-> z1
  amat["z1","y"] <- 2; amat["y","z1"] <- 3; # z1 -> y
  amat["z2","y"] <- 2; amat["y","z2"] <- 3; # z2 -> y
  amat["z2","w"] <- 1; amat["w","z2"] <- 1; # w o-o z2
  #plotAG(amat)

  x <- c("x")
  y <- c("y")
  z <- c("z1", "z2")
  return(list(amat=amat, x=x, y=y, z=z))
}


##############################################################################
# Asia / Lung Cancer (hypothetical model)                                    #
# Source: bnlearn                                                            #
# Paper: https://www.eecis.udel.edu/~shatkay/Course/papers/Lauritzen1988.pdf #
##############################################################################

# A DAG represented by the PAG in Fig. 2c) in R50 - https://causalai.net/r50.pdf
getTrueAsiaDAG <- function(lat = NULL, ret_dagg = TRUE, ret_bnfit = TRUE) {
  # TRUE DAG
  allvars <-  c("dysp", "either", "smoke", "bronc", "tub", "lung", "asia", "xray")
  p <- length(allvars)
  true.dag <- matrix(0, p, p)
  colnames(true.dag) <- rownames(true.dag) <- allvars
  true.dag["asia", "tub"] <- 0; true.dag["tub", "asia"] <- 1; # asia -> tub
  true.dag["tub", "either"] <- 0; true.dag["either", "tub"] <- 1; # tub -> either
  true.dag["either", "xray"] <- 0; true.dag["xray", "either"] <- 1; # either -> xray
  true.dag["either", "dysp"] <- 0; true.dag["dysp", "either"] <- 1; # either -> dysp
  true.dag["lung", "either"] <- 0; true.dag["either", "lung"] <- 1; # lung -> either
  true.dag["smoke", "lung"] <- 0; true.dag["lung", "smoke"] <- 1; # smoke -> lung
  true.dag["smoke", "bronc"] <- 0; true.dag["bronc", "smoke"] <- 1; # smoke -> bronc
  true.dag["bronc", "dysp"] <- 0; true.dag["dysp", "bronc"] <- 1; # bronc -> dysp
  # plot(as(t(true.dag), "graphNEL"))

  x <- c("either")
  y <- c("dysp")
  z <- c("smoke") #"bronc"
  lat <- lat

  dagg=NULL
  if (ret_dagg) {
    dagg <- pcalg::pcalg2dagitty(true.dag, colnames(true.dag),type="dag")
    #plot(dagg)

    dagitty::latents(dagg) <- lat
    dagitty::exposures(dagg) <- x
    dagitty::outcomes(dagg) <- y

    edges <- dagitty::edges(dagg)
    edges[which(edges$v == "asia" & edges$w=="tub"), "beta"] <- 0.15 # runif(1, 0.2, 0.3)
    edges[which(edges$v == "tub" & edges$w=="either"), "beta"] <- 0.5 # runif(1, 0.6, 0.8)
    edges[which(edges$v == "lung" & edges$w=="either"), "beta"] <- 0.5 # runif(1, 0.6, 0.8)
    edges[which(edges$v == "smoke" & edges$w=="bronc"), "beta"] <- 0.7 # runif(1, 0.6, 0.8)
    edges[which(edges$v == "smoke" & edges$w=="lung"), "beta"] <- 0.7 # runif(1, 0.6, 0.8)
    edges[which(edges$v == "either" & edges$w=="xray"), "beta"] <- 0.4 # runif(1, 0.6, 0.8)
    edges[which(edges$v == "bronc" & edges$w=="dysp"), "beta"] <- 0.3 # runif(1, 0.6, 0.8)
    edges[which(edges$v == "either" & edges$w=="dysp"), "beta"] <- 0.55 # runif(1, 0.6, 0.8)

    dagg_str <- as.character(dagg)
    for (eid in 1:nrow(edges)) {
      edge_str <- paste0(edges[eid, "v"], " ", edges[eid, "e"], " ", edges[eid, "w"])
      new_edge_str <- paste0(edge_str, " [beta = ", edges[eid, "beta"], "]")
      dagg_str <- gsub(edge_str, new_edge_str, dagg_str)
    }
    dagg <- dagitty::dagitty(dagg_str)


    # R <- impliedCovarianceMatrix(dagg, b.default = NULL, b.lower = -0.6,
    #                              b.upper = 0.6, eps = 1, standardized = TRUE)
    # R[R < 1e-16] <- 0
    # R < 0
    # R <- round(R, 7)
    # R < 0.1 & R > 0
    #
    # library(matrixcalc)
    # is.positive.definite(R, tol=1e-16)

  }

  #plot(getPAG(dagg))

  if (ret_bnfit) {
    asia_fit <- getCustomAsiaBN() # readRDS("asia.rds")
    # Note
    # true.dag <- t(amat(asia_fit)) # it needs to be transposed to meet pcalg notation
    # dagg <-  pcalg2dagitty(true.dag, colnames(true.dag),type="dag")
    # plot(as(t(true.dat), "graphNEL")) # DAG
    # graphviz.chart(asia_fit)
  }


  return(list(amat=true.dag, x=x, y=y, z=z, lat=lat, dagg=dagg, bnfit=asia_fit))
}

getTrueAsiaPAG <- function() {
  true.dag <- getTrueAsiaDAG(lat = NULL, ret_dagg=TRUE)
  true_fci_out <- getPAG(true.dag$dagg)
  true.pag <- as(true_fci_out, "amat")
  temp <- matrix(true.pag, ncol(true.dag$amat), ncol(true.dag$amat))
  colnames(temp) <- rownames(temp) <- colnames(true.pag)
  true.pag <- temp

  #plotAG(true.pag)

  x <- true.dag$x
  y <- true.dag$y
  z <- true.dag$z
  return(list(amat=true.pag, x=x, y=y, z=z))
}

