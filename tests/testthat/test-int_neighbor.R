context("int_neighbor")

#load data
colkas <- qtl::read.cross(format="csvs",dir="./",genfile="ColKas_geno.csv",phefile = "ColKas_pheno.csv",
                          na.strings = c("_"), estimate.map=TRUE, crosstype = "riself")
colkas <- colkas[1:2,1:30]
colkas_genoprob <- qtl::calc.genoprob(colkas, step=4)
x <- colkas$pheno[,2]
y <- colkas$pheno[,3]
smap_colkas <- data.frame(x,y)

#F2
set.seed(1234)
data("fake.f2",package="qtl")
fake_f2 <- fake.f2[1:2,1:30]
smap_f2 <- cbind(runif(qtl::nind(fake_f2),1,100),runif(qtl::nind(fake_f2),1,100))
genoprobs_f2 <- qtl::calc.genoprob(fake_f2,step=4)

#backcross
set.seed(1234)
data("fake.bc",package="qtl")
fake_bc <- fake.bc[1:2,1:30]
smap_bc <- cbind(runif(qtl::nind(fake_bc),1,100),runif(qtl::nind(fake_bc),1,100))
genoprobs_bc <- qtl::calc.genoprob(fake_bc,step=4)

test_that(
  desc = "LOD_range",
  code = {
    colkas_int <- int_neighbor(genoprobs=colkas_genoprob,
                               pheno=log(colkas$pheno[,5]+1),
                               contrasts=c(TRUE,FALSE,TRUE),
                               smap=smap_colkas, scale=7,
                               addcovar=as.matrix(colkas$pheno[,7:9]),
                               addQTL=c("c1_nga280"), intQTL="c1_nga280")

    f2_int <- int_neighbor(genoprobs=genoprobs_f2,
                           pheno=fake_f2$pheno[,1],
                           contrasts=c(TRUE,TRUE,TRUE),
                           smap=smap_f2, scale=20,
                           addcovar=as.matrix(fake_f2$pheno$sex),
                           addQTL=c("c1_D1M318","c1_D1M212"), intQTL="c1_D1M212",
                           grouping=fake_f2$pheno$sex)

    bc_int <- int_neighbor(genoprobs=genoprobs_bc,
                           pheno=fake_bc$pheno[,1],
                           contrasts=c(TRUE,TRUE,FALSE),
                           smap=smap_bc, scale=59,
                           addcovar=as.matrix(cbind(fake_bc$pheno$sex,fake_bc$pheno$age)),
                           addQTL=c("c1_D1M318","c1_D1M212"), intQTL="c1_D1M212")


    expect_true(all(colkas_int$LOD_self>=0))
    expect_true(all(f2_int$LOD_self>=0))
    expect_true(all(bc_int$LOD_self>=0))

    expect_true(all(colkas_int$LOD_nei>=0))
    expect_true(all(f2_int$LOD_nei>=0))
    expect_true(all(bc_int$LOD_nei>=0))

  }
)
