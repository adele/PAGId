library(pcalg)
library(dagitty)
library(igraph)
library(causaleffect)
library(bnlearn)
library(PAGId)
source("PAGInstances.R")
source("ProbUtils.R")
source("EstimateCausalEffect.R")

#######################################
# Creating PAG of Fig. 5 in the paper #
#######################################

adag <- getDAG2n()
apag <- getPAG(adag$dagg) # generates the PAG using the oracle
#plot(apag) # PAG
amat <- apag@amat
x <- adag$x
y <- adag$y
z <- adag$z

plot(as(t(adag$amat), "graphNEL")) # DAG showing even the latent variables
dev.off()
plot(apag) # PAG
dev.off()
plotAG(amat) # PAG from the adj matrix

#################################
# Identifying the causal effect #
#################################

retPAG <- CIDP(amat, x, y, z)
retPAG$id
print(paste0("ID from PAG: ", retPAG$Qexpr[[retPAG$query]]))


#############################################################
# Estimating the causal effect from data and the ID formula #
#############################################################


N=10000
dat <- generateDataset(adag$dagg, N, type="binary", alpha_cutoff=0.1, verbose=FALSE)

QopElem = retPAG$Qop[[retPAG$query]]
Qop = retPAG$Qop

vnames <- c(adag$x, adag$z)
vals_grid <- expand.grid(comprehenr::to_list(for (var in vnames) c(0,1)))
colnames(vals_grid) <- vnames
vals_grid[[adag$y]] <- 1

other_vnames <- setdiff(colnames(amat), colnames(vals_grid))
other_vals <- matrix(0, nrow=nrow(vals_grid), ncol=length(other_vnames))
colnames(other_vals) <- other_vnames
vals_grid <- cbind(vals_grid, other_vals)


py.dox.z_PAG <- c()
for (i in 1:nrow(vals_grid)) {
  curvals <- unlist(vals_grid[i,])
  py.dox.z_PAG <- c(py.dox.z_PAG, computeQop(dat, QopElem, Qop, values=curvals)$prob)
}
py.dox.z_PAG <- cbind(vals_grid, "py_dox.z"=py.dox.z_PAG)



