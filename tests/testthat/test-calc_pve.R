context("calc_pve")

#load data
colkas <- qtl::read.cross(format="csvs",dir="./",genfile="ColKas_geno.csv",phefile = "ColKas_pheno.csv",
                          na.strings = c("_"), estimate.map=TRUE, crosstype = "riself")
colkas <- qtl2::convert2cross2(colkas)
gmap_colkas <- qtl2::insert_pseudomarkers(colkas$gmap, step=2)
colkas_genoprob <- qtl2::calc_genoprob(colkas,gmap_colkas)
x <- colkas$pheno[,2]
y <- colkas$pheno[,3]
smap_colkas <- data.frame(x,y)
s_colkas <- quantile(dist(smap_colkas),c(0.1*(0:10)))

#F2
set.seed(1234)
data("fake.f2",package="qtl")
fake_f2 <- qtl2::convert2cross2(fake.f2)
fake_f2 <- subset(fake_f2,chr=c(1:19))
smap_f2 <- cbind(runif(qtl2::n_ind(fake_f2),1,100),runif(qtl2::n_ind(fake_f2),1,100))
gmap_f2 <- qtl2::insert_pseudomarkers(fake_f2$gmap, step=2)
genoprobs_f2 <- qtl2::calc_genoprob(fake_f2,gmap_f2)
s_f2 <- quantile(dist(smap_f2),c(0.1*(1:10)))

#backcross
set.seed(1234)
data("fake.bc",package="qtl")
fake_bc <- qtl2::convert2cross2(fake.bc)
fake_bc <- subset(fake_bc,chr=c(1:19))
smap_bc <- cbind(runif(qtl2::n_ind(fake_bc),1,100),runif(qtl2::n_ind(fake_bc),1,100))
gmap_bc <- qtl2::insert_pseudomarkers(fake_bc$gmap, step=2)
genoprobs_bc <- qtl2::calc_genoprob(fake_bc,gmap_bc)
s_bc <- quantile(dist(smap_bc),c(0.1*(1:10)))

f2_bin <- as.numeric(fake_f2$pheno[,1]>mean(fake_f2$pheno[,1]))
bc_bin <- as.numeric(fake_bc$pheno[,1]>mean(fake_bc$pheno[,1]))

test_that(
  desc = "pve_range",
  code = {

    colkas_pve <- calc_pve(genoprobs=colkas_genoprob,
                           pheno=log(colkas$pheno[,4]+1),
                           gmap=gmap_colkas, contrasts=c(TRUE,FALSE,TRUE),
                           smap=smap_colkas, s_seq=s_colkas,
                           addcovar=colkas$pheno[,6:8], fig=FALSE
                           )

    f2_pve <- calc_pve(genoprobs=genoprobs_f2,
                       pheno=fake_f2$pheno[,1],
                       gmap=gmap_f2, contrasts=c(TRUE,TRUE,TRUE),
                       smap=smap_f2, s_seq=s_f2,
                       addcovar=as.matrix(fake_f2$covar), fig=FALSE
                       )

    bc_pve <- calc_pve(genoprobs=genoprobs_bc,
                       pheno=fake_bc$pheno[,1],
                       gmap=gmap_bc, contrasts=c(TRUE,TRUE,FALSE),
                       smap=smap_bc, s_seq=s_bc,
                       addcovar=as.matrix(fake_bc$covar), fig=FALSE
                       )

    expect_true(all(round(colkas_pve[,2],1)>=0))
    expect_true(all(round(f2_pve[,2],1)>=0))
    expect_true(all(round(bc_pve[,2],1)>=0))

    expect_true(all(round(colkas_pve[,2],1)<=1))
    expect_true(all(round(f2_pve[,2],1)<=1))
    expect_true(all(round(bc_pve[,2],1)<=1))

})

test_that(
  desc = "binary_action",
  code = {

    colkas_pveBin <- calc_pve(genoprobs=colkas_genoprob,
                              gmap=gmap_colkas, contrasts=c(TRUE,FALSE,TRUE),
                              pheno=colkas$pheno[,5],
                              smap=smap_colkas,s_seq=s_colkas,
                              addcovar=NULL, response="binary", fig=FALSE
                              )

    f2_pveBin <- calc_pve(genoprobs=genoprobs_f2,
                          pheno=f2_bin,
                          gmap=gmap_f2, contrasts=c(TRUE,TRUE,TRUE),
                          smap=smap_f2, s_seq=s_f2,
                          addcovar=as.matrix(fake_f2$covar), response="binary", fig=FALSE
                          )

    bc_pveBin <- calc_pve(genoprobs=genoprobs_bc,
                          pheno=bc_bin,
                          gmap=gmap_bc, contrasts=c(TRUE,TRUE,FALSE),
                          smap=smap_bc, s_seq=s_bc,
                          addcovar=as.matrix(fake_bc$covar), response="binary", fig=FALSE
                          )

    expect_true(all(is.na(colkas_pveBin[,3])))
    expect_true(all(is.na(f2_pveBin[,3])))
    expect_true(all(is.na(bc_pveBin[,3])))

  }
)
