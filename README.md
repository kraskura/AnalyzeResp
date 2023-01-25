
<!-- README.md is generated from README.Rmd. Please edit that file -->

# AnalyzeResp

<!-- badges: start -->

[![R-CMD-check](https://github.com/kraskura/AnalyzeResp/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/kraskura/AnalyzeResp/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

This project is in transition to update. The older version is documented
here: <https://kraskura.github.io/AnalyzeResp_0/> Please note, the older
version will not be changed or improved.

## Installation:

The development version of AnalyzeResp can be installed from
[GitHub](https://github.com/) with:

``` r
install.packages("devtools") # if devtools is not installed
devtools::install_github("kraskura/AnalyzeResp")
library(AnalyzeResp)
```

## What has changed since the last version:

- The folder names and process
- Multiple SMR and SDA files combine on the spot
- Convert units in the convert text function
- The ‘csv_files’ folder is ‘outside’
- SMR inventory lives ‘outside’ in the main wd (as set at the start)
- Dont’ need to reset wd
- Names of the functions

## To do list:

- Clean the code and get syntax under control. (some arguments have “.”
  in names, while some have “\_” in their names). Carryover of old not
  good practices of coding
- Make sure the imported libraries load
- Rework the manual and instructions
- Rework an example code with the following scenarios
  - One files and multiple (glued files)
  - Various background strategies (one back prior, only back pre, etc)
- Document the independent sliding window analysis
- Use this for Loligo respos (Witrox file recorded files)
- Revise exported data files and make sure each column of data is
  necessary, revise accordingly.
- Improve and label plots
- Revise and delete the out-commented unused code
- Use lintr to clean up old mess
