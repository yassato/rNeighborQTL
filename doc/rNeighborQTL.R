## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,  fig.width = 4, fig.height = 4,
  comment = "#>"
)

## ----input--------------------------------------------------------------------
colkas <- qtl::read.cross(format="csvs",dir="../inst",
                          genfile="ColKas_geno.csv",
                          phefile = "ColKas_pheno.csv",
                          na.strings = c("_"), estimate.map=TRUE, crosstype = "riself"
                          )

colkas_genoprob <- qtl::calc.genoprob(colkas, step=2)

## ----pve----------------------------------------------------------------------
library(rNeighborQTL)
x <- colkas$pheno[,2]
y <- colkas$pheno[,3]
smap_colkas <- data.frame(x,y)

s_seq <- quantile(dist(smap_colkas),c(0.1*(1:10)))
colkas_pve <- calc_pve(genoprobs=colkas_genoprob,
                       pheno=log(colkas$pheno[,5]+1),
                       smap=smap_colkas, s_seq=s_seq,
                       addcovar=as.matrix(colkas$pheno[,7:9]) 
                       )

## ----eff, fig.width=4, fig.height=8-------------------------------------------
colkas_eff <- eff_neighbor(genoprobs=colkas_genoprob,
                           pheno=log(colkas$pheno[,5]+1),
                           smap=smap_colkas, scale=7,
                           addcovar=as.matrix(colkas$pheno[,7:9])
                           )

## ----LOD----------------------------------------------------------------------
colkas_scan <- scan_neighbor(genoprobs=colkas_genoprob, 
                             pheno=log(colkas$pheno[,5]+1),
                             smap=smap_colkas, scale=7, 
                             addcovar=as.matrix(colkas$pheno[,7:9])
                             )
plot_nei(colkas_scan)

## ----perm, eval=FALSE---------------------------------------------------------
#  colkas_perm <- perm_neighbor(genoprobs=colkas_genoprob,
#                               pheno=log(colkas$pheno[,5]+1),
#                               smap=smap_colkas, scale=7,
#                               addcovar=as.matrix(colkas$pheno[,6:8]),
#                               times=3, p_val=c(0.5,0.1)
#                               )

## ----self---------------------------------------------------------------------
plot_nei(colkas_scan, type="self")
colkas_scanone <- qtl::scanone(colkas_genoprob,
                            pheno.col=log(colkas$pheno$holes+1),
                            addcovar=as.matrix(colkas$pheno[,7:9]),
                            method="hk")
plot(colkas_scanone)

## ----CIM, eval=FALSE----------------------------------------------------------
#  colkas_cim <- scan_neighbor(genoprobs=colkas_genoprob,
#                              pheno=log(colkas$pheno[,5]+1),
#                              smap=smap_colkas, scale=7,
#                              addcovar=as.matrix(colkas$pheno[,7:9]),
#                              addQTL="c4_nga8"
#                              )
#  plot_nei(colkas_cim)

## ----int, eval=FALSE----------------------------------------------------------
#  colkas_int <- int_neighbor(genoprobs=colkas_genoprob,
#                             pheno=log(colkas$pheno[,5]+1),
#                             smap=smap_colkas, scale=7,
#                             addcovar=as.matrix(colkas$pheno[,7:9]),
#                             addQTL="c4_nga8", intQTL="c4_nga8"
#                             )
#  
#  plot_nei(colkas_int, type="int")

## ----bin----------------------------------------------------------------------
s_seq <- quantile(dist(smap_colkas),c(0.1*(1:10)))
colkas_pveBin <- calc_pve(genoprobs=colkas_genoprob, 
                          pheno=colkas$pheno[,7],
                          smap=smap_colkas, s_seq=s_seq,
                          response="binary", addcovar=as.matrix(colkas$pheno[,8:9]), 
                          fig=TRUE)

colkas_scanBin <- scan_neighbor(genoprobs=colkas_genoprob, 
                                pheno=colkas$pheno[,7],
                                smap=smap_colkas, scale=2.24,
                                addcovar=as.matrix(colkas$pheno[,8:9]), 
                                response="binary")

plot_nei(colkas_scanBin)

## ----fake---------------------------------------------------------------------
#F2 lines
set.seed(1234)
data("fake.f2",package="qtl")
fake_f2 <- subset(fake.f2, chr=1:19)
smap_f2 <- cbind(runif(qtl::nind(fake_f2),1,100),runif(qtl::nind(fake_f2),1,100))
genoprobs_f2 <- qtl::calc.genoprob(fake_f2,step=2)
s_seq <- quantile(dist(smap_f2),c(0.1*(1:10)))

nei_eff <- sim_nei_qtl(genoprobs_f2, a2=0.5, d2=0.5, 
                       smap=smap_f2, 
                       scale=s_seq[1], n_QTL=1
                       )

pve_f2 <- calc_pve(genoprobs=genoprobs_f2,
                   pheno=nei_eff$nei_y,
                   smap=smap_f2, s_seq=s_seq[1:5],
                   addcovar=as.matrix(cbind(fake_f2$pheno$sex,fake_f2$pheno$pgm)),
                   fig=FALSE)
    
deltaPVE <- pve_f2[-1,3] - c(0,pve_f2[1:4,3])
argmax_s <- s_seq[1:5][deltaPVE==max(deltaPVE)]
    
scan_f2 <- scan_neighbor(genoprobs=genoprobs_f2,
                         pheno=nei_eff$nei_y,
                         smap=smap_f2, scale=argmax_s,
                         addcovar=as.matrix(cbind(fake_f2$pheno$sex,fake_f2$pheno$pgm))
                         )
    
plot_nei(scan_f2)

## ----bc-----------------------------------------------------------------------
#backcross lines
set.seed(1234)
data("fake.bc",package="qtl")
fake_bc <- subset(fake.bc, chr=1:19)
smap_bc <- cbind(runif(qtl::nind(fake_bc),1,100),runif(qtl::nind(fake_bc),1,100))
genoprobs_bc <- qtl::calc.genoprob(fake_bc,step=2)
s_seq <- quantile(dist(smap_bc),c(0.1*(1:10)))

nei_eff <- sim_nei_qtl(genoprobs_bc, a2=0.3, d2=-0.3, 
                       smap=smap_bc, 
                       scale=s_seq[1], n_QTL=1)

pve_bc <- calc_pve(genoprobs=genoprobs_bc,
                   pheno=nei_eff$nei_y,
                   smap=smap_bc, s_seq=s_seq[1:5],
                   addcovar=as.matrix(cbind(fake_bc$pheno$sex,fake_bc$pheno$age)),
                   fig=FALSE)
    
deltaPVE <- pve_bc[-1,3] - c(0,pve_bc[1:4,3])
argmax_s <- s_seq[1:5][deltaPVE==max(deltaPVE)]
    
scan_bc <- scan_neighbor(genoprobs=genoprobs_bc,
                         pheno=nei_eff$nei_y,
                         smap=smap_bc, scale=argmax_s,
                         addcovar=as.matrix(cbind(fake_bc$pheno$sex,fake_bc$pheno$age))
                         )

plot_nei(scan_bc)

