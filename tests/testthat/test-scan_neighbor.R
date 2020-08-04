context("scan_neighbor")

#load data
colkas <- qtl::read.cross(format="csvs",dir="./",genfile="ColKas_geno.csv",phefile = "ColKas_pheno.csv",
                          na.strings = c("_"), estimate.map=TRUE, crosstype = "riself")
colkas <- colkas[1:2,1:30]
colkas_genoprob <- qtl::calc.genoprob(colkas, step=4)
x <- colkas$pheno[,2]
y <- colkas$pheno[,3]
smap_colkas <- data.frame(x,y)
s_colkas <- quantile(dist(smap_colkas),c(0.1*(0:10)))

#F2
set.seed(1234)
data("fake.f2",package="qtl")
fake_f2 <- fake.f2[1:2,1:30]
smap_f2 <- cbind(runif(qtl::nind(fake_f2),1,100),runif(qtl::nind(fake_f2),1,100))
genoprobs_f2 <- qtl::calc.genoprob(fake_f2,step=4)
s_f2 <- quantile(dist(smap_f2),c(0.1*(1:10)))

#backcross
set.seed(1234)
data("fake.bc",package="qtl")
fake_bc <- fake.bc[1:2,1:30]
smap_bc <- cbind(runif(qtl::nind(fake_bc),1,100),runif(qtl::nind(fake_bc),1,100))
genoprobs_bc <- qtl::calc.genoprob(fake_bc,step=4)
s_bc <- quantile(dist(smap_bc),c(0.1*(1:10)))

f2_bin <- as.numeric(fake_f2$pheno[,1]>mean(fake_f2$pheno[,1]))
bc_bin <- as.numeric(fake_bc$pheno[,1]>mean(fake_bc$pheno[,1]))

test_that(
  desc="LOD_range",
  code = {
    colkas_scan <- scan_neighbor(genoprobs=colkas_genoprob,
                                 pheno=log(colkas$pheno[,5]+1),
                                 contrasts=c(TRUE,FALSE,TRUE),
                                 smap=smap_colkas, scale=7,
                                 addcovar=as.matrix(colkas$pheno[,7:9])
                                 )

    f2_scan <- scan_neighbor(genoprobs=genoprobs_f2,
                             pheno=fake_f2$pheno[,1],
                             contrasts=c(TRUE,TRUE,TRUE),
                             smap=smap_f2, scale=19.37,
                             addcovar=as.matrix(fake_f2$pheno$sex)
                             )

    bc_scan <- scan_neighbor(genoprobs=genoprobs_bc,
                             pheno=fake_bc$pheno[,1],
                             contrasts=c(TRUE,TRUE,FALSE),
                             smap=smap_bc, scale=59,
                             addcovar=as.matrix(cbind(fake_bc$pheno$sex,fake_bc$pheno$age))
                             )

    expect_true(all(colkas_scan$LOD_self>=0))
    expect_true(all(f2_scan$LOD_self>=0))
    expect_true(all(bc_scan$LOD_self>=0))

    expect_true(all(colkas_scan$LOD_nei>=0))
    expect_true(all(f2_scan$LOD_nei>=0))
    expect_true(all(bc_scan$LOD_nei>=0))
  }
)

test_that(
  desc = "self_equal",
  code = {
    colkas_scan <- scan_neighbor(genoprobs=colkas_genoprob,
                                 pheno=log(colkas$pheno[,5]+1),
                                 contrasts=c(TRUE,FALSE,TRUE),
                                 smap=smap_colkas, scale=7,
                                 addcovar=as.matrix(colkas$pheno[,7:9])
                                 )

    f2_scan <- scan_neighbor(genoprobs=genoprobs_f2,
                             pheno=fake_f2$pheno[,1],
                             contrasts = c(TRUE,TRUE,TRUE),
                             smap=smap_f2, scale=19.37,
                             addcovar=as.matrix(fake_f2$pheno$sex)
                             )

    bc_scan <- scan_neighbor(genoprobs=genoprobs_bc,
                             pheno=fake_bc$pheno[,1],
                             contrasts=c(TRUE,TRUE,FALSE),
                             smap=smap_bc, scale=59,
                             addcovar=as.matrix(cbind(fake_bc$pheno$sex,fake_bc$pheno$age))
                             )

    colkas_scan1 <- qtl::scanone(colkas_genoprob,pheno.col=log(colkas$pheno[,5]+1),addcovar=as.matrix(colkas$pheno[,7:9]),method="hk")
    f2_scan1 <- qtl::scanone(genoprobs_f2,pheno.col=fake_f2$pheno[,1],addcovar=as.matrix(fake_f2$pheno$sex),method="ehk")
    bc_scan1 <- qtl::scanone(genoprobs_bc,pheno.col=fake_bc$pheno[,1],addcovar=as.matrix(cbind(fake_bc$pheno$sex,fake_bc$pheno$age)),method="hk")

    expect_equal(round(stats::cor(colkas_scan$LOD_self, colkas_scan1[1:nrow(colkas_scan),3]),1),1)
    expect_equal(round(stats::cor(f2_scan$LOD_self, f2_scan1[1:nrow(f2_scan),3]),1),1)
    expect_equal(round(stats::cor(bc_scan$LOD_self, bc_scan1[1:nrow(bc_scan),3]),1),1)
  }
)

test_that(
  desc = "CIM_equal",
  code = {
    colkas_scan2 <- scan_neighbor(genoprobs=colkas_genoprob,
                                  pheno=log(colkas$pheno[,5]+1),
                                  contrasts=c(TRUE,FALSE,TRUE),
                                  smap=smap_colkas, scale=7,
                                  addcovar=as.matrix(colkas$pheno[,7:9]),
                                  addQTL="c1_nga280")

    f2_scan2 <- scan_neighbor(genoprobs=genoprobs_f2,
                              pheno=fake_f2$pheno[,1],
                              contrasts = c(TRUE,TRUE,TRUE),
                              smap=smap_f2, scale=19.37,
                              addcovar=as.matrix(fake_f2$pheno$sex),
                              addQTL=c("c1_D1M318","c1_D1M212"))

    bc_scan2 <- scan_neighbor(genoprobs=genoprobs_bc,
                              pheno=fake_bc$pheno[,1],
                              contrasts=c(TRUE,TRUE,FALSE),
                              smap=smap_bc, scale=59,
                              addcovar=as.matrix(cbind(fake_bc$pheno$sex,fake_bc$pheno$age)),
                              addQTL=c("c1_D1M318","c1_D1M212"))

    expect_equal(colkas_scan2["c1_nga280",4],0)
    expect_equal(f2_scan2["c1_D1M318",4],0)
    expect_equal(f2_scan2["c1_D1M212",4],0)
    expect_equal(bc_scan2["c1_D1M318",4],0)
    expect_equal(bc_scan2["c1_D1M212",4],0)
  }
)

test_that(
  desc = "bin_LOD_range",
  code = {
    colkas_bin <- scan_neighbor(genoprobs=colkas_genoprob,
                                pheno=colkas$pheno[,6],
                                contrasts=c(TRUE,FALSE,TRUE),
                                smap=smap_colkas, scale=7,
                                addcovar=as.matrix(colkas$pheno[,7:9]),
                                response="binary")

    f2_bin <- scan_neighbor(genoprobs=genoprobs_f2,
                            pheno=f2_bin,
                            contrasts=c(TRUE,TRUE,TRUE),
                            smap=smap_f2, scale=19.37,
                            addcovar=as.matrix(fake_f2$pheno$sex),
                            response="binary")

    bc_bin <- scan_neighbor(genoprobs=genoprobs_bc,
                            pheno=bc_bin,
                            contrasts=c(TRUE,TRUE,FALSE),
                            smap=smap_bc, scale=59,
                            addcovar=as.matrix(cbind(fake_bc$pheno$sex,fake_bc$pheno$age)),
                            response="binary")

    expect_true(all(colkas_bin$LOD_self>=0))
    expect_true(all(f2_bin$LOD_self>=0))
    expect_true(all(bc_bin$LOD_self>=0))

    expect_true(all(colkas_bin$LOD_nei>=0))
    expect_true(all(f2_bin$LOD_nei>=0))
    expect_true(all(bc_bin$LOD_nei>=0))
  }
)
