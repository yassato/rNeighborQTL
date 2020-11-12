# rNeighborQTL: An R package for interval QTL mapping of neighbor effects.  
The "rNeighborQTL" package provides a set of functions to perform a genome scan of neighbor QTL effects.  

## Installation
CRAN version is available at https://cran.r-project.org/package=rNeighborQTL.  
To install the deveper version via GitHub, run *devtools::install_github("yassato/rNeighborQTL", repo="master")*.  

## Dependency
Note that the rNeighborQTL requires the following R packages.  
- gaston
- parallel
- qtl

## Release Note  
version 1.1.1 (developer version): partial PVEs separated in calc_PVEnei(); plot_nei() refined.  
version 1.1.0 (CRAN release version): calc_neiprob(), eff_neighbor(), scan_neighbor() refactored.  
version 1.0.0: Initial version registered in CRAN.  
