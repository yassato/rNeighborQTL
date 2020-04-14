context("eff_neighbor")

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

test_that(
  desc = "skip_AB_riself",
  code = {

    colkas_eff1 <- eff_neighbor(genoprobs=colkas_genoprob,
                                pheno=log(colkas$pheno[,4]+1),
                                gmap=gmap_colkas, contrasts=c(TRUE,FALSE,TRUE),
                                smap=smap_colkas, scale=7, addcovar=colkas$pheno[,6:8], fig=FALSE)
    expect_true(all(is.na(colkas_eff1$d1)))
    expect_true(all(is.na(colkas_eff1$d2)))

  }
)

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
  desc = "symmetry_bc",
  code = {
    effect_bc <- eff_neighbor(genoprobs=genoprobs_bc,
                              pheno=fake_bc$pheno[,1],
                              gmap=gmap_bc, contrasts=c(TRUE,TRUE,FALSE),
                              smap=smap_bc, scale=59, addcovar=as.matrix(fake_bc$covar), fig=FALSE)
    expect_equal(effect_bc$a1, effect_bc$d1)
  }
)
