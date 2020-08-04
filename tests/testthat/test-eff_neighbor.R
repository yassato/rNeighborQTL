context("eff_neighbor")

#load data
colkas <- qtl::read.cross(format="csvs",dir="./",genfile="ColKas_geno.csv",phefile = "ColKas_pheno.csv",
                          na.strings = c("_"), estimate.map=TRUE, crosstype = "riself")
colkas <- colkas[1:2,1:50]
colkas_genoprob <- qtl::calc.genoprob(colkas, step=4)
x <- colkas$pheno[,2]
y <- colkas$pheno[,3]
smap_colkas <- data.frame(x,y)
s_colkas <- quantile(dist(smap_colkas),c(0.1*(0:10)))

test_that(
  desc = "skip_AB_riself",
  code = {
    colkas_eff1 <- eff_neighbor(genoprobs=colkas_genoprob,
                                pheno=log(colkas$pheno[,5]+1),
                                smap=smap_colkas, scale=7,
                                addcovar=as.matrix(colkas$pheno[,7:9]), fig=FALSE)
    expect_true(all(is.na(colkas_eff1$d1)))
    expect_true(all(is.na(colkas_eff1$d2)))

  }
)

#backcross
data("fake.bc",package="qtl")
fake_bc <- fake.bc[1:2,1:50]
smap_bc <- cbind(runif(qtl::nind(fake_bc),1,100),runif(qtl::nind(fake_bc),1,100))
genoprobs_bc <- qtl::calc.genoprob(fake_bc,step=4)

test_that(
  desc = "symmetry_bc",
  code = {
    effect_bc <- eff_neighbor(genoprobs=genoprobs_bc,
                              pheno=fake_bc$pheno[,1],
                              smap=smap_bc, scale=59,
                              addcovar=as.matrix(cbind(fake_bc$pheno$sex,fake_bc$pheno$age)),
                              fig=FALSE)
    expect_equal(effect_bc$a1, effect_bc$d1)
  }
)
