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

defaultPath <- getwd()
stopifnot(length(grep("RUnit", defaultPath)) == 1)

temp <- list.files(file.path(defaultPath, "..", "..", "pkg", "dirichlet", "R"), full.names=TRUE)
junk <- lapply(temp, source)

testSuite.d <-
  defineTestSuite("dirichlet", dirs=file.path(defaultPath),
                  testFileRegexp="^runit_[[:alnum:]]+.[rR]$")

testResult <- runTestSuite(testSuite.d)

htmlFile <- file.path(defaultPath, "Test Results.html")

## warning expected about gcc compiler
suppressWarnings(printHTMLProtocol(testResult, fileName=htmlFile))

browseURL(htmlFile, browser=getOption("browser"))

