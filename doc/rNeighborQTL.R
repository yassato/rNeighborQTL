## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,  fig.width = 4, fig.height = 4,
  comment = "#>"
)

## ----input--------------------------------------------------------------------
colkas <- qtl::read.cross(format="csvs",dir="../inst",
                    genfile="ColKas_geno.csv",
                    phefile = "ColKas_pheno.csv",
                    na.strings = c("_"), estimate.map=TRUE, crosstype = "riself")
colkas <- qtl2::convert2cross2(colkas)

gmap_colkas <- qtl2::insert_pseudomarkers(colkas$gmap, step=2)
colkas_genoprob <- qtl2::calc_genoprob(colkas,gmap_colkas)

## ----pve----------------------------------------------------------------------
library(rNeighborQTL)
x <- colkas$pheno[,2]
y <- colkas$pheno[,3]
smap_colkas <- data.frame(x,y)

s_seq <- quantile(dist(smap_colkas),c(0.1*(1:10)))
colkas_pve <- calc_pve(genoprobs=colkas_genoprob,
                      pheno=log(colkas$pheno[,4]+1),
                      gmap=gmap_colkas, contrasts=c(TRUE,FALSE,TRUE),
                      addcovar=colkas$pheno[,6:8], 
                      smap=smap_colkas, s_seq=s_seq
                      )

## ----eff, fig.width=4, fig.height=8-------------------------------------------
colkas_eff <- eff_neighbor(genoprobs=colkas_genoprob,
                           pheno=log(colkas$pheno[,4]+1),
                           gmap=gmap_colkas, contrasts=c(TRUE,FALSE,TRUE),
                           smap=smap_colkas, scale=7,
                           addcovar=colkas$pheno[,6:8]
                           )

## ----LOD----------------------------------------------------------------------
colkas_scan <- scan_neighbor(genoprobs=colkas_genoprob, 
                             pheno=log(colkas$pheno[,4]+1),
                             gmap=gmap_colkas, contrasts=c(TRUE,FALSE,TRUE),
                             smap=smap_colkas, scale=7, 
                             addcovar=colkas$pheno[,6:8]
                             )
plot_nei(colkas_scan)

## ----perm---------------------------------------------------------------------
colkas_perm <- perm_neighbor(genoprobs=colkas_genoprob, pheno=log(colkas$pheno[,4]+1),
                            gmap=gmap_colkas, contrasts=c(TRUE,FALSE,TRUE),
                            smap=smap_colkas, scale=7,
                            addcovar=colkas$pheno[,6:8],
                            times=99, p_val=c(0.1,0.05,0.01))
print(colkas_perm)

## ----self---------------------------------------------------------------------
plot_nei(colkas_scan, type="self")
colkas_scan1 <- qtl2::scan1(colkas_genoprob,pheno=colkas$pheno[,2])
plot(colkas_scan1, map=gmap_colkas)

## ----CIM----------------------------------------------------------------------
colkas_cim <- scan_neighbor(genoprobs=colkas_genoprob, pheno=log(colkas$pheno[,4]+1),
                            gmap=gmap_colkas, contrasts=c(TRUE,FALSE,TRUE),
                            smap=smap_colkas, scale=7,
                            addcovar=colkas$pheno[,6:8],
                            addQTL="nga8"
                            )
plot_nei(colkas_cim)

## ----int----------------------------------------------------------------------
colkas_int <- int_neighbor(genoprobs=colkas_genoprob, 
                           pheno=log(colkas$pheno[,4]+1), 
                           gmap=gmap_colkas, contrasts=c(TRUE,FALSE,TRUE),
                           smap=smap_colkas, scale=7, 
                           addcovar=colkas$pheno[,6:8], 
                           addQTL="nga8", intQTL="nga8"
                           )
plot_nei(colkas_int, type="int")

## ----bin, fig.width=4, fig.height=8-------------------------------------------
s_seq <- quantile(dist(smap_colkas),c(0.1*(1:10)))
colkas_pveBin <- calc_pve(genoprobs=colkas_genoprob, pheno=colkas$pheno[,5],
                       contrasts=c(TRUE,FALSE,TRUE), addcovar=NULL,
                       gmap=gmap_colkas, smap=smap_colkas,
                       s_seq=s_seq, response="binary", fig=FALSE
                       )

colkas_scanBin <- scan_neighbor(genoprobs=colkas_genoprob, pheno=colkas$pheno[,5],
                                gmap=gmap_colkas, contrasts=c(TRUE,FALSE,TRUE),
                                smap_colkas, scale=7.82, response="binary")

par(mfcol=c(2,1))
plot_nei(colkas_scanBin,type="self")
plot_nei(colkas_scanBin,type="neighbor")

## ----fake---------------------------------------------------------------------
#demo using F2
set.seed(1234)
data("fake.f2",package="qtl")
fake_f2 <- qtl2::convert2cross2(fake.f2)
fake_f2 <- subset(fake_f2,chr=c(1:19))
smap_f2 <- cbind(runif(qtl2::n_ind(fake_f2),1,100),runif(qtl2::n_ind(fake_f2),1,100))
gmap_f2 <- qtl2::insert_pseudomarkers(fake_f2$gmap, step=2)
genoprobs_f2 <- qtl2::calc_genoprob(fake_f2,gmap_f2)
s_seq <- quantile(dist(smap_f2),c(0.1*(1:10)))
pve_f2 <- calc_pve(genoprobs=genoprobs_f2,
                     pheno=fake_f2$pheno[,1], gmap=gmap_f2,
                     smap=smap_f2, s_seq=s_seq,
                     contrasts=c(TRUE,TRUE,TRUE),
                     addcovar=as.matrix(fake_f2$covar), fig=FALSE
                     )
scan_f2 <- scan_neighbor(genoprobs=genoprobs_f2,
                         pheno=fake_f2$pheno[,1], gmap=gmap_f2,
                         contrasts = c(TRUE,TRUE,TRUE), smap=smap_f2,
                         scale=19.37, 
                         addcovar=as.matrix(fake_f2$covar)
                         )
plot_nei(scan_f2)

## ----bc-----------------------------------------------------------------------
#backcross lines
set.seed(1234)
data("fake.bc",package="qtl")
fake_bc <- qtl2::convert2cross2(fake.bc)
fake_bc <- subset(fake_bc,chr=c(1:19))
smap_bc <- cbind(runif(qtl2::n_ind(fake_bc),1,100),runif(qtl2::n_ind(fake_bc),1,100))
s_seq <- quantile(dist(smap_bc),c(0.1*(1:10)))
gmap_bc <- qtl2::insert_pseudomarkers(fake_bc$gmap, step=2)
genoprobs_bc <- qtl2::calc_genoprob(fake_bc,gmap_bc)
pve_bc <- calc_pve(genoprobs=genoprobs_bc,pheno=fake_bc$pheno[,1],
                   gmap=gmap_bc, smap=smap_bc, s_seq=s_seq,
                   contrasts=c(TRUE,TRUE,FALSE),
                   addcovar=as.matrix(fake_bc$covar), fig=FALSE
                   )
scan_bc <- scan_neighbor(genoprobs=genoprobs_bc,
                        pheno=fake_bc$pheno[,1],
                        gmap=gmap_bc, contrasts = c(TRUE,TRUE,FALSE),
                        smap=smap_bc, scale=59,
                        addcovar=as.matrix(fake_bc$covar))
plot_nei(scan_bc)

