context("calc_pve")

#load data
colkas <- qtl::read.cross(format="csvs",dir="./",genfile="ColKas_geno.csv",phefile = "ColKas_pheno.csv",
                          na.strings = c("_"), estimate.map=TRUE, crosstype = "riself")
colkas_genoprob <- qtl::calc.genoprob(colkas, step=2)
x <- colkas$pheno[,2]
y <- colkas$pheno[,3]
smap_colkas <- data.frame(x,y)
s_colkas <- quantile(dist(smap_colkas),c(0.1*(0:10)))

#F2
set.seed(1234)
data("fake.f2",package="qtl")
fake_f2 <- fake.f2[1:19,]
smap_f2 <- cbind(runif(qtl::nind(fake_f2),1,100),runif(qtl::nind(fake_f2),1,100))
genoprobs_f2 <- qtl::calc.genoprob(fake_f2,step=2)
s_f2 <- quantile(dist(smap_f2),c(0.1*(1:10)))

#backcross
set.seed(1234)
data("fake.bc",package="qtl")
fake_bc <- fake.bc[1:19,]
smap_bc <- cbind(runif(qtl::nind(fake_bc),1,100),runif(qtl::nind(fake_bc),1,100))
genoprobs_bc <- qtl::calc.genoprob(fake_bc,step=2)
s_bc <- quantile(dist(smap_bc),c(0.1*(1:10)))

f2_bin <- as.numeric(fake_f2$pheno[,1]>mean(fake_f2$pheno[,1]))
bc_bin <- as.numeric(fake_bc$pheno[,1]>mean(fake_bc$pheno[,1]))

test_that(
  desc = "pve_range",
  code = {

    colkas_pve <- calc_pve(genoprobs=colkas_genoprob,
                           pheno=log(colkas$pheno[,5]+1),
                           contrasts=c(TRUE,FALSE,TRUE),
                           smap=smap_colkas, s_seq=s_colkas[1:5],
                           addcovar=as.matrix(colkas$pheno[,7:9]),
                           fig=FALSE)

    f2_pve <- calc_pve(genoprobs=genoprobs_f2,
                       pheno=fake_f2$pheno[,1],
                       contrasts=c(TRUE,TRUE,TRUE),
                       smap=smap_f2, s_seq=s_f2[1:5],
                       addcovar=as.matrix(cbind(fake_f2$pheno$sex,fake_f2$pheno$pgm)),
                       fig=FALSE)

    bc_pve <- calc_pve(genoprobs=genoprobs_bc,
                       pheno=fake_bc$pheno[,1],
                       contrasts=c(TRUE,TRUE,FALSE),
                       smap=smap_bc, s_seq=s_bc[1:5],
                       addcovar=as.matrix(cbind(fake_bc$pheno$sex,fake_bc$pheno$pgm)),
                       fig=FALSE)

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
                              contrasts=c(TRUE,FALSE,TRUE),
                              pheno=colkas$pheno[,6],
                              smap=smap_colkas,s_seq=s_colkas[1:5],
                              addcovar=NULL,
                              response="binary", fig=FALSE)

    f2_pveBin <- calc_pve(genoprobs=genoprobs_f2,
                          pheno=f2_bin,
                          contrasts=c(TRUE,TRUE,TRUE),
                          smap=smap_f2, s_seq=s_f2[1:5],
                          addcovar=as.matrix(cbind(fake_f2$pheno$sex,fake_f2$pheno$pgm)),
                          response="binary", fig=FALSE)

    bc_pveBin <- calc_pve(genoprobs=genoprobs_bc,
                          pheno=bc_bin,
                          contrasts=c(TRUE,TRUE,FALSE),
                          smap=smap_bc, s_seq=s_bc[1:5],
                          addcovar=as.matrix(cbind(fake_bc$pheno$sex,fake_bc$pheno$pgm)),
                          response="binary", fig=FALSE)

    expect_true(all(is.na(colkas_pveBin[,3])))
    expect_true(all(is.na(f2_pveBin[,3])))
    expect_true(all(is.na(bc_pveBin[,3])))

  }
)
