context("perm_neighbor")

#load data
colkas <- qtl::read.cross(format="csvs",dir="./",genfile="ColKas_geno.csv",phefile = "ColKas_pheno.csv",
                          na.strings = c("_"), estimate.map=TRUE, crosstype = "riself")
colkas <- colkas[1:2,1:50]
colkas_genoprob <- qtl::calc.genoprob(colkas, step=4)
x <- colkas$pheno[,2]
y <- colkas$pheno[,3]
smap_colkas <- data.frame(x,y)
s_colkas <- quantile(dist(smap_colkas),c(0.1*(0:10)))

#F2
set.seed(1234)
data("fake.f2",package="qtl")
fake_f2 <- fake.f2[1:2,1:50]
smap_f2 <- cbind(runif(qtl::nind(fake_f2),1,100),runif(qtl::nind(fake_f2),1,100))
genoprobs_f2 <- qtl::calc.genoprob(fake_f2,step=4)
s_f2 <- quantile(dist(smap_f2),c(0.1*(1:10)))

#backcross
set.seed(1234)
data("fake.bc",package="qtl")
fake_bc <- fake.bc[1:2,1:50]
smap_bc <- cbind(runif(qtl::nind(fake_bc),1,100),runif(qtl::nind(fake_bc),1,100))
genoprobs_bc <- qtl::calc.genoprob(fake_bc,step=4)
s_bc <- quantile(dist(smap_bc),c(0.1*(1:10)))

test_that(
  desc = "perm_LOD_min",
  code = {
    perm_colkas <- perm_neighbor(genoprobs=colkas_genoprob,
                                 pheno=log(colkas$pheno[,5]+1),
                                 smap=smap_colkas, scale=7,
                                 addcovar=as.matrix(colkas$pheno[,7:9]),
                                 times=2, p_val=1.0, type="neighbor")

    perm_f2 = perm_neighbor(genoprobs=genoprobs_f2,
                            pheno=fake_f2$pheno[,1],
                            smap=smap_f2, scale=28.2,
                            addcovar=as.matrix(fake_f2$pheno$sex),
                            times=2, p_val=1.0, type="neighbor")

    perm_bc = perm_neighbor(genoprobs=genoprobs_bc,
                            pheno=fake_bc$pheno[,1],
                            smap=smap_bc, scale=50.7,
                            addcovar=as.matrix(cbind(fake_bc$pheno$sex,fake_bc$pheno$age)),
                            times=2, p_val=1.0, type="neighbor")

    expect_true(perm_colkas>=0)
    expect_true(perm_f2>=0)
    expect_true(perm_bc>=0)
  }
)

test_that(
  desc = "int_perm_LOD_min",
  code = {
    perm_colkas_int <- perm_neighbor(genoprobs=colkas_genoprob,
                               pheno=log(colkas$pheno[,5]+1),
                               smap=smap_colkas, scale=7,
                               addcovar=as.matrix(colkas$pheno[,7:9]),
                               addQTL=c("c1_nga280"), intQTL="c1_nga280",
                               times=2, p_val=1.0, type="int")

    perm_f2_int <- perm_neighbor(genoprobs=genoprobs_f2,
                           pheno=fake_f2$pheno[,1],
                           smap=smap_f2, scale=20,
                           addcovar=as.matrix(fake_f2$pheno$sex),
                           addQTL=c("c1_D1M318","c1_D1M212"), intQTL="c1_D1M212",
                           grouping=fake_f2$pheno$sex,
                           times=2, p_val=1.0, type="int")

    perm_bc_int <- perm_neighbor(genoprobs=genoprobs_bc,
                           pheno=fake_bc$pheno[,1],
                           smap=smap_bc, scale=59,
                           addcovar=as.matrix(cbind(fake_bc$pheno$sex,fake_bc$pheno$age)),
                           addQTL=c("c1_D1M318","c1_D1M212"), intQTL="c1_D1M212",
                           times=2, p_val=1.0, type="int")

    expect_true(perm_colkas_int>=0)
    expect_true(perm_f2_int>=0)
    expect_true(perm_bc_int>=0)
  }
)
