
setwd("C:\\Users\\User\\Documents\\ZueF\\Hainbuche_oktoploid\\code_information\\oktoploid_analysis\\polyfreq_test")

remove.packages("ZueFPolyPloid")

install.packages("C:\\Users\\User\\Documents\\ZueF\\Hainbuche_oktoploid\\code_information\\oktoploid_analysis\\polyfreq", repos = NULL, type = "source")
library(polyfreq)
library(polysat)
library(testthat)


####Testcase for mygen Dataset, with 2 Loci and Ploidie 4 ##

test_that("de_silva_freq_EM matches polysat::deSilvaFreq on mygen (self=0.2)", {
  
  load("mygen.Rdata")

  #check the calculations with the new R package
  start <- Sys.time()
  calc_freq_polyfreq <- polyfreq::de_silva_freq_EM(mygen, self = 0.2)
  end <- Sys.time()
  polyfreq_timediff <- end - start  # Returns difftime object

  ####
  #check the calculations with the new R package
  start <- Sys.time()
  calc_freq_polysat <- polysat::deSilvaFreq(mygen, self = 0.2)
  end <- Sys.time()
  polysat_timediff <- end - start  # Returns difftime object

  # Numeric frequencies must be identical to polysat within machine precision
  expect_equal(calc_freq_polyfreq, calc_freq_polysat, tolerance = 1e-10)

})

############################################################

#### Testcase for Tetr Dataset, with 10 Loci, up to 15 distinct alleles per Locus
####  and Ploidie 4 ##


test_that("de_silva_freq_EM matches polysat::deSilvaFreq on Tetr (self=0.5)", {

  load("Tetr.Rdata")

  Ploidies(Tetr) <- 4

  #check the calculations with the new R package
  start <- Sys.time()
  calc_freq_polyfreq <- polyfreq::de_silva_freq_EM(Tetr, self = 0.5)
  end <- Sys.time()
  polyfreq_timediff <- end - start  # Returns difftime object

  ####
  #check the calculations with the new R package
  start <- Sys.time()
  calc_freq_polysat <- polysat::deSilvaFreq(Tetr, self = 0.5)
  end <- Sys.time()
  polysat_timediff <- end - start  # Returns difftime object
 
  expect_equal(calc_freq_polyfreq, calc_freq_polysat, tolerance = 1e-10)

})
############################################################


############################################################

#### Testcase for Hornbeam Dataset, with 8 Loci, up to 25 distinct alleles per Locus
####  and Ploidie 8 ##
#This does not work at the moment for both calculation approaches, as
#the matrices get too big
#better ways of coping with the algorithm are needed

test_that("de_silva_freq_EM matches polysat::deSilvaFreq on hornbeam_genambig (self=0.5)", {
  skip("This does not work at the moment for both calculation approaches, as the matrices get too big")
  
  load("hornbeam_genambig.Rdata")

  #check the calculations with the new R package
  start <- Sys.time()
  calc_freq_polyfreq <- polyfreq::de_silva_freq_EM(hornbeam_genambig, self = 0.5)
  end <- Sys.time()
  polyfreq_timediff <- end - start  # Returns difftime object

  ####
  #check the calculations with the new R package
  start <- Sys.time()
  calc_freq_polysat <- polysat::deSilvaFreq(hornbeam_genambig, self = 0.5)
  end <- Sys.time()
  polysat_timediff <-end - start  # Returns difftime object

  expect_equal(calc_freq_polyfreq, calc_freq_polysat, tolerance = 1e-10)

})

##########################################################