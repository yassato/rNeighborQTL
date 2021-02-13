# rNeighborQTL  
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
version 1.1.2 (developer version): rebuild using R version 4.0.3; citations updated.  
version 1.1.1 (CRAN release version): partial PVEs separated in calc_PVEnei(); plot_nei() refined.  
version 1.1.0: calc_neiprob(), eff_neighbor(), scan_neighbor() refactored.  
version 1.0.0: Initial version registered in CRAN.  

## References  
Sato Y, Takeda K, Nagano AJ (2021) Neighbor QTL: an interval mapping method for quantitative trait loci underlying plant neighborhood effects. G3; Genes|Genomes|Genetics https://doi.org/10.1093/g3journal/jkab017  
