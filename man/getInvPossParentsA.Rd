% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PAGUtils.R
\name{getInvPossParentsA}
\alias{getInvPossParentsA}
\title{Gets the "invisible" possible parents of A, i.e.,
the set \{ V :  V o-> A or V -> A with an invisible edge\}}
\usage{
getInvPossParentsA(amat_v, amat, nodeA)
}
\arguments{
\item{amat_v}{TODO}

\item{amat}{the \code{pcalg}-compatible adjacency matrix of the PAG}

\item{nodeA}{name of the variable A}
}
\value{
names of all "invisible" possible parents of A
}
\description{
Gets the "invisible" possible parents of A, i.e.,
the set \{ V :  V o-> A or V -> A with an invisible edge\}
}
