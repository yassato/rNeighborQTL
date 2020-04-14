context("scan_neighbor")

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
  desc="LOD_range",
  code = {
    colkas_scan <- scan_neighbor(genoprobs=colkas_genoprob,
                                 pheno=log(colkas$pheno[,4]+1),
                                 gmap=gmap_colkas, contrasts=c(TRUE,FALSE,TRUE),
                                 smap=smap_colkas, scale=7,
                                 addcovar=colkas$pheno[,6:8])

    f2_scan <- scan_neighbor(genoprobs=genoprobs_f2,
                             pheno=fake_f2$pheno[,1],
                             gmap=gmap_f2, contrasts = c(TRUE,TRUE,TRUE),
                             smap=smap_f2, scale=19.37,
                             addcovar=as.matrix(fake_f2$covar))

    bc_scan <- scan_neighbor(genoprobs=genoprobs_bc,
                             pheno=fake_bc$pheno[,1],
                             gmap=gmap_bc, contrasts=c(TRUE,TRUE,FALSE),
                             smap=smap_bc, scale=59,
                             addcovar=as.matrix(fake_bc$covar))

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
                                 pheno=log(colkas$pheno[,4]+1),
                                 gmap=gmap_colkas, contrasts=c(TRUE,FALSE,TRUE),
                                 smap=smap_colkas, scale=7,
                                 addcovar=colkas$pheno[,6:8])

    f2_scan <- scan_neighbor(genoprobs=genoprobs_f2,
                             pheno=fake_f2$pheno[,1],
                             gmap=gmap_f2, contrasts = c(TRUE,TRUE,TRUE),
                             smap=smap_f2, scale=19.37,
                             addcovar=as.matrix(fake_f2$covar))

    bc_scan <- scan_neighbor(genoprobs=genoprobs_bc,
                             pheno=fake_bc$pheno[,1],
                             gmap=gmap_bc, contrasts=c(TRUE,TRUE,FALSE),
                             smap=smap_bc, scale=59,
                             addcovar=as.matrix(fake_bc$covar))

    colkas_scan1 <- qtl2::scan1(colkas_genoprob,pheno=log(colkas$pheno[,4]+1),addcovar=colkas$pheno[,6:8])
    f2_scan1 <- qtl2::scan1(genoprobs_f2,pheno=fake_f2$pheno[,1],addcovar=fake_f2$covar)
    bc_scan1 <- qtl2::scan1(genoprobs_bc,pheno=fake_bc$pheno[,1],addcovar=fake_bc$covar)

    expect_equal(round(stats::cor(colkas_scan$LOD_self, colkas_scan1[1:nrow(colkas_scan),]),2),1)
    expect_equal(round(stats::cor(f2_scan$LOD_self, f2_scan1[1:nrow(f2_scan),]),2),1)
    expect_equal(round(stats::cor(bc_scan$LOD_self, bc_scan1[1:nrow(bc_scan),]),2),1)
  }
)

test_that(
  desc = "CIM_equal",
  code = {
    colkas_scan2 <- scan_neighbor(genoprobs=colkas_genoprob,
                                  pheno=log(colkas$pheno[,4]+1),
                                  gmap=gmap_colkas, contrasts=c(TRUE,FALSE,TRUE),
                                  smap=smap_colkas, scale=7, addcovar=colkas$pheno[,6:8],
                                  addQTL="nga8")

    f2_scan2 <- scan_neighbor(genoprobs=genoprobs_f2,
                             pheno=fake_f2$pheno[,1],
                             gmap=gmap_f2, contrasts = c(TRUE,TRUE,TRUE),
                             smap=smap_f2, scale=19.37,
                             addcovar=as.matrix(fake_f2$covar),
                             addQTL=c("D1M318","D1M212"))

    bc_scan2 <- scan_neighbor(genoprobs=genoprobs_bc,
                             pheno=fake_bc$pheno[,1],
                             gmap=gmap_bc, contrasts=c(TRUE,TRUE,FALSE),
                             smap=smap_bc, scale=59,
                             addcovar=as.matrix(fake_bc$covar),
                             addQTL=c("D1M318","D1M212"))

    expect_equal(colkas_scan2["nga8",4],0)
    expect_equal(f2_scan2["D1M318",4],0)
    expect_equal(f2_scan2["D1M212",4],0)
    expect_equal(bc_scan2["D1M318",4],0)
    expect_equal(bc_scan2["D1M212",4],0)
  }
)

test_that(
  desc = "bin_LOD_range",
  code = {
    colkas_bin <- scan_neighbor(genoprobs=colkas_genoprob,
                                pheno=colkas$pheno[,5],
                                gmap=gmap_colkas, contrasts=c(TRUE,FALSE,TRUE),
                                smap=smap_colkas, scale=7,
                                addcovar=colkas$pheno[,6:8], response="binary")

    f2_bin <- scan_neighbor(genoprobs=genoprobs_f2,
                            pheno=f2_bin,
                            gmap=gmap_f2, contrasts = c(TRUE,TRUE,TRUE),
                            smap=smap_f2, scale=19.37,
                            addcovar=as.matrix(fake_f2$covar), response="binary")

    bc_bin <- scan_neighbor(genoprobs=genoprobs_bc,
                            pheno=bc_bin,
                            gmap=gmap_bc, contrasts=c(TRUE,TRUE,FALSE),
                            smap=smap_bc, scale=59,
                            addcovar=as.matrix(fake_bc$covar), response="binary")

    expect_true(all(colkas_bin$LOD_self>=0))
    expect_true(all(f2_bin$LOD_self>=0))
    expect_true(all(bc_bin$LOD_self>=0))

    expect_true(all(colkas_bin$LOD_nei>=0))
    expect_true(all(f2_bin$LOD_nei>=0))
    expect_true(all(bc_bin$LOD_nei>=0))
  }
)
