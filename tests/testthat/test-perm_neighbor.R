context("perm_neighbor")

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

test_that(
  desc = "perm_LOD_min",
  code = {
    perm_colkas <- perm_neighbor(genoprobs=colkas_genoprob,
                                 pheno=log(colkas$pheno[,4]+1),
                                 gmap=gmap_colkas, contrasts=c(TRUE,FALSE,TRUE),
                                 smap=smap_colkas, scale=7, addcovar=colkas$pheno[,6:8],
                                 times=9, p_val=1.0, type="neighbor")

    perm_f2 = perm_neighbor(genoprobs=genoprobs_f2,
                            pheno=fake_f2$pheno[,1],
                            gmap=gmap_f2, contrasts = c(TRUE,TRUE,TRUE),
                            smap=smap_f2, scale=28.2, addcovar=as.matrix(fake_f2$covar),
                            times=9, p_val=1.0, type="neighbor")

    perm_bc = perm_neighbor(genoprobs=genoprobs_bc,
                            pheno=fake_bc$pheno[,1],
                            gmap=gmap_bc, contrasts = c(TRUE,TRUE,FALSE),
                            smap=smap_bc, scale=50.7, addcovar=as.matrix(fake_bc$covar),
                            times=9, p_val=1.0, type="neighbor")

    expect_true(perm_colkas>=0)
    expect_true(perm_f2>=0)
    expect_true(perm_bc>=0)
  }
)

test_that(
  desc = "int_perm_LOD_min",
  code = {
    perm_colkas_int <- perm_neighbor(genoprobs=colkas_genoprob,
                               pheno=log(colkas$pheno[,4]+1),
                               gmap=gmap_colkas, contrasts=c(TRUE,FALSE,TRUE),
                               smap=smap_colkas, scale=7, addcovar=colkas$pheno[,6:8],
                               addQTL=c("nga8"), intQTL="nga8",
                               times=9, p_val=1.0, type="int")

    perm_f2_int <- perm_neighbor(genoprobs=genoprobs_f2,
                           pheno=fake_f2$pheno[,1],
                           gmap=gmap_f2, contrasts=c(TRUE,TRUE,TRUE),
                           smap=smap_f2, scale=20,
                           addcovar=as.matrix(fake_f2$covar),
                           addQTL=c("D1M318","D1M212"), intQTL="D1M212",
                           grouping=fake_f2$covar$pgm,
                           times=9, p_val=1.0, type="int")

    perm_bc_int <- perm_neighbor(genoprobs=genoprobs_bc,
                           pheno=fake_bc$pheno[,1],
                           gmap=gmap_bc, contrasts = c(TRUE,TRUE,FALSE),
                           smap=smap_bc, scale=59,
                           addcovar=NULL,
                           addQTL=c("D1M318","D1M212"), intQTL="D1M212",
                           times=9, p_val=1.0, type="int")

    expect_true(perm_colkas_int>=0)
    expect_true(perm_f2_int>=0)
    expect_true(perm_bc_int>=0)
  }
)
