context("plot_nei")

#F2
set.seed(1234)
data("fake.f2",package="qtl")
fake_f2 <- qtl2::convert2cross2(fake.f2)
fake_f2 <- subset(fake_f2,chr=c(1:19))
smap_f2 <- cbind(runif(qtl2::n_ind(fake_f2),1,100),runif(qtl2::n_ind(fake_f2),1,100))
gmap_f2 <- qtl2::insert_pseudomarkers(fake_f2$gmap, step=2)
genoprobs_f2 <- qtl2::calc_genoprob(fake_f2,gmap_f2)
s_f2 <- quantile(dist(smap_f2),c(0.1*(1:10)))

test_that(
  desc = "plot_error_catch",
  code = {
    f2_scan <- scan_neighbor(genoprobs=genoprobs_f2,
                             pheno=fake_f2$pheno[,1],
                             gmap=gmap_f2, contrasts = c(TRUE,TRUE,TRUE),
                             smap=smap_f2, scale=19.37,
                             addcovar=as.matrix(fake_f2$covar))

    f2_int <- int_neighbor(genoprobs=genoprobs_f2,
                           pheno=fake_f2$pheno[,1],
                           gmap=gmap_f2, contrasts=c(TRUE,TRUE,TRUE),
                           smap=smap_f2, scale=20,
                           addcovar=as.matrix(fake_f2$covar),
                           addQTL=c("D1M318","D1M212"), intQTL="D1M212",
                           grouping=fake_f2$covar$pgm)

    expect_error(plot_nei(f2_scan,type="int"))
    expect_error(plot_nei(f2_int,type="self"))
    expect_error(plot_nei(f2_int,type="neighbor"))

  }
)
