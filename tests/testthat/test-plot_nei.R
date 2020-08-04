context("plot_nei")

#F2
set.seed(1234)
data("fake.f2",package="qtl")
fake_f2 <- fake.f2[1:2,1:30]
smap_f2 <- cbind(runif(qtl::nind(fake_f2),1,100),runif(qtl::nind(fake_f2),1,100))
genoprobs_f2 <- qtl::calc.genoprob(fake_f2,step=4)


test_that(
  desc = "plot_error_catch",
  code = {
    f2_scan <- scan_neighbor(genoprobs=genoprobs_f2,
                             pheno=fake_f2$pheno[,1],
                             contrasts = c(TRUE,TRUE,TRUE),
                             smap=smap_f2, scale=19.37,
                             addcovar=as.matrix(fake_f2$pheno$sex)
                             )

    f2_int <- int_neighbor(genoprobs=genoprobs_f2,
                           pheno=fake_f2$pheno[,1],
                           contrasts=c(TRUE,TRUE,TRUE),
                           smap=smap_f2, scale=20,
                           addcovar=as.matrix(fake_f2$pheno$sex),
                           addQTL=c("c1_D1M318","c1_D1M212"), intQTL="c1_D1M212",
                           grouping=fake_f2$pheno$sex)

    expect_error(plot_nei(f2_scan,type="int"))
    expect_error(plot_nei(f2_int,type="self"))
    expect_error(plot_nei(f2_int,type="neighbor"))

  }
)
