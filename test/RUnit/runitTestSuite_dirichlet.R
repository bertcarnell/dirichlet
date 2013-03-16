################################################################################
#
# Program Name:  runitTestSuite_dirichlet.R
# Purpose:       To provide test functions for the dirichlet package
# Author:        Rob Carnell
# Date:          July 2007
#
# Required Packages:  RUnit
# R version:          2.4.0 (>=2.3.1)
#
################################################################################

rm(list=ls())
require(RUnit)

## used as part of package
#require(dirichlet)
#defaultPath <- chartr("/", "//", paste(.path.package("dirichlet"), "/RUnit", sep=""))

################# used in development ##########################################
defaultPath <- file.path("c:", "Documents and Settings", "carnellr",
                         "My Documents", "Repositories", "BTRA",
                         "Trunk", "Source", "Rscripts", "dirichlet")
                         
temp <- list.files(file.path(defaultPath, "R"), full.names=TRUE)
junk <- lapply(temp, source)

testSuite.d <-
  defineTestSuite("dirichlet", dirs=file.path(defaultPath, "RUnit"),
                  testFileRegexp="^runit_[[:alnum:]]+.[rR]$")

testResult <- runTestSuite(testSuite.d)

################# used in development ##########################################

htmlFile <- file.path(defaultPath, "RUnit", "Test Results.html")

## warning expected about gcc compiler
suppressWarnings(printHTMLProtocol(testResult, fileName=htmlFile))

browseURL(htmlFile, browser=getOption("browser"))

