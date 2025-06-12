# PAGId R Package: Algorithms for (Conditional) Causal Effect Identification in Partial Ancestral Graphs

## Overview

This package implements the CIDP and IDP algorithms for identifying (conditional) causal effects from a Partial Ancentral Graph (PAG). Technical details are provided in the NeurIPS 2022 paper by Jaber A., Ribeiro A. H., Zhang J., and Bareinboim E., (2022) entitled ["Causal Identification under Markov equivalence: Calculus, Algorithm, and Completeness"](https://proceedings.neurips.cc/paper_files/paper/2022/hash/17a9ab4190289f0e1504bbb98d1d111a-Abstract-Conference.html).

## Requirements:

First, install R (>= 3.5.0) and the following necessary R packages:
```r
install.packages(c("dagitty", "pcalg"), dependencies=TRUE)
```

If you also want to run the simulations, then the following R packages are also necessary:
```r
install.packages(c("bnlearn", "causaleffect", "igraph"), dependencies=TRUE)
```

## Installation:

We can now proceed with the installation of the PAGId R package. 


You can download the latest tar.gz file with the source code of the PAGId R package, available at https://github.com/adele/PAGId/releases/latest, and install it with the following command, where path_to_file represents the full path and file name of the tar.gz file:

```r
install.packages(path_to_file, repos=NULL, type="source", dependencies=TRUE)
```

Or you can install the development version directly from GitHub. Make sure you have the devtools R package installed. If not, install it with install.packages("devtools", dependencies=TRUE).

```r
devtools::install_github("adele/PAGId", dependencies=TRUE)
```
