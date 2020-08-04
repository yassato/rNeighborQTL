context("decompose_genoprobs")

#load data
colkas <- qtl::read.cross(format="csvs",dir="./",
                          genfile="ColKas_geno.csv",phefile = "ColKas_pheno.csv",
                          na.strings = c("_"), estimate.map=TRUE, crosstype = "riself")
colkas <- colkas[1:2,1:50]
colkas_genoprob <- qtl::calc.genoprob(colkas, step=4)

#F2
set.seed(1234)
data("fake.f2",package="qtl")
fake_f2 <- fake.f2[1:2,1:50]
genoprobs_f2 <- qtl::calc.genoprob(fake_f2,step=4)

#backcross
set.seed(1234)
data("fake.bc",package="qtl")
fake_bc <- fake.bc[1:2,1:50]
genoprobs_bc <- qtl::calc.genoprob(fake_bc,step=4)

test_that(
  desc = "prob_sum_1",
  code = {
    geno_colkas <- decompose_genoprobs(colkas_genoprob, contrasts=NULL)
    geno_f2 <- decompose_genoprobs(genoprobs_f2, contrasts=NULL)
    geno_bc <- decompose_genoprobs(genoprobs_bc, contrasts=NULL)

    expect_equal(sum(round(geno_colkas$AA+geno_colkas$BB,3)!=1), 0)
    expect_equal(sum(round(geno_f2$AA+geno_f2$AB+geno_f2$BB,3)!=1), 0)
    expect_equal(sum(round(geno_bc$AA+geno_bc$AB,3)!=1), 0)
})
