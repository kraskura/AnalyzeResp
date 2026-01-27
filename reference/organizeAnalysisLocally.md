# convenient function to create a set of local working directories to streamline workflow

Creates a collection of local directories that help organizing the
analysis by exporting output files to specific folders associated with
each metabolic performance function (not required for smooth data
analysis).

## Usage

``` r
organizeAnalysisLocally(
  SMR.folder = TRUE,
  MMR.folder = TRUE,
  SDA.folder = TRUE,
  MMR_SMR_AS_EPOC.folder = TRUE,
  BACK_RESP.folder = TRUE
)
```

## Arguments

- SMR.folder:

  create a local SMR folder, default = TRUE

- MMR.folder:

  create a local MMR folder, default = TRUE

- SDA.folder:

  create a local SDA folder, default = TRUE

- MMR_SMR_AS_EPOC.folder:

  create a local MMR_SMR_AS_EPOC.folder, default = TRUE

- BACK_RESP.folder:

  create a local BACK_RESP.folder, default = TRUE

## Value

The output from [`print`](https://rdrr.io/r/base/print.html)
