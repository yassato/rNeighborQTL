## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----input_files--------------------------------------------------------------
colkas <- qtl::read.cross(format="csvs",dir=".",
                    genfile="../inst/ColKasGeno2001.csv",
                    phefile = "../inst/ColKas_fakePheno.csv",
                    na.strings = c("_"), estimate.map=TRUE)

colkas <- qtl2::convert2cross2(colkas)

## ----get_selfprob-------------------------------------------------------------
gmap_colkas <- qtl2::insert_pseudomarkers(colkas$gmap, step=1)
colkas_genoprob <- qtl2::calc_genoprob(colkas,gmap_colkas,error_prob=0.002)

## ----neiprob------------------------------------------------------------------
library(rNeighborQTL)
set.seed(1)
x <- runif(qtl2::n_ind(colkas),1,100)
y <- runif(qtl2::n_ind(colkas),1,100)
smap_colkas <- data.frame(x,y)

colkas_neiprob <- calc_neiprob(genoprobs=colkas_genoprob,
                              gmap=gmap_colkas, contrasts=c(TRUE,TRUE,TRUE),
                              smap=smap_colkas, scale=33,
                              a2=1, d2=0)

## ----calc_pve-----------------------------------------------------------------
s_seq <- quantile(dist(smap_colkas),c(0.1*(1:10))) #get quantiles for pairwise distances
colkas_pve <- calc_pve(genoprobs=colkas_genoprob, pheno=colkas$pheno[,2],
                       gmap=gmap_colkas, contrasts=c(TRUE,TRUE,TRUE),
                       smap=smap_colkas, s_seq=s_seq)

## ----eff_neighbor, fig.width=4, fig.height=8----------------------------------
colkas_effect1 <- eff_neighbor(genoprobs=colkas_genoprob,
                                  pheno=colkas$pheno[,2],
                                  gmap=gmap_colkas, contrasts=c(TRUE,TRUE,TRUE),
                                  smap=smap_colkas, scale=33, fig=TRUE)

## ----scan_neighbor------------------------------------------------------------
colkas_scan1 <- scan_neighbor(genoprobs=colkas_genoprob,
                              pheno=colkas$pheno[,2],
                              gmap=gmap_colkas, contrasts=c(TRUE,TRUE,TRUE),
                              smap=smap_colkas, scale=33)
plot_nei(colkas_scan1, type="neighbor")

## ----perm_neighbor------------------------------------------------------------
perm_colkas <- perm_neighbor(genoprobs=colkas_genoprob, pheno=colkas$pheno[,2],
                              gmap=gmap_colkas, contrasts=c(TRUE,TRUE,TRUE),
                              smap=smap_colkas, scale=33,
                              times=9, p_val=0.01)
plot_nei(colkas_scan1, type="neighbor", th=perm_colkas$LOD_th)

## ----bonf_neighbor------------------------------------------------------------
plot_nei(colkas_scan1, type="neighbor", th="bonf",ylim=c(0,4.5))

## ----self---------------------------------------------------------------------
plot_nei(colkas_scan1, type="self")
colkas_scan1 <- qtl2::scan1(colkas_genoprob,pheno=colkas$pheno[,2])
plot(colkas_scan1, map=gmap_colkas)

## ----CIM----------------------------------------------------------------------
colkas_scan2 <- scan_neighbor(genoprobs=colkas_genoprob,
                              pheno=colkas$pheno[,2],
                              gmap=gmap_colkas, contrasts=c(TRUE,TRUE,TRUE),
                              smap=smap_colkas, scale=33,
                              addQTL=c("nga59","nga280"))
plot_nei(colkas_scan2)

## ----glm----------------------------------------------------------------------
colkas_scan3 <- scan_neighbor(genoprobs=colkas_genoprob,
                              pheno=colkas$pheno[,3],
                              gmap=gmap_colkas, contrasts=c(TRUE,TRUE,TRUE),
                              smap=smap_colkas, scale=33,
                              response="count")
plot_nei(colkas_scan3)

## ----inbred-------------------------------------------------------------------
colkas <- qtl::read.cross(format="csvs",dir=".",
                      genfile="../inst/ColKasGeno2001_inbred.csv",
                      phefile = "../inst/ColKas_fakePheno.csv",
                      na.strings = c("_"), estimate.map=TRUE,
                      crosstype = "riself")

colkas <- qtl2::convert2cross2(colkas)
gmap_colkas <- qtl2::insert_pseudomarkers(colkas$gmap, step=1)
colkas_genoprob <- qtl2::calc_genoprob(colkas,gmap_colkas,error_prob=0.002)

colkas_scan4 <- scan_neighbor(genoprobs=colkas_genoprob,
                              pheno=colkas$pheno[,3],
                              gmap=gmap_colkas, contrasts=c(TRUE,FALSE,TRUE),
                              smap=smap_colkas, scale=33,
                              response="count")
plot_nei(colkas_scan4)

## ----fake---------------------------------------------------------------------
?scan_neighbor

