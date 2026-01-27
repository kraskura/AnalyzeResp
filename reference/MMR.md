# Estimates regressions for MMR.

Analyzes O~2~ and temperature data recorded to obtain regressions
representative of individual's maximum metabolic rate. The function can
handle 4 measurement cycles.

## Usage

``` r
MMR(
  data.MMR,
  cycles,
  cycle_start,
  cycle_end,
  mmr_Ch1 = 0,
  mmr_Ch2 = 0,
  mmr_Ch3 = 0,
  mmr_Ch4 = 0,
  clean_Ch1 = c(0, 0),
  clean_Ch2 = c(0, 0),
  clean_Ch3 = c(0, 0),
  clean_Ch4 = c(0, 0),
  N_Ch = 4,
  local_path = TRUE,
  date_format = c("%m/%d/%Y %H:%M:%S", "GMT"),
  inventory_data = NULL
)
```

## Arguments

- data.MMR:

  The name of the raw input file; .csv, a character.

- cycles:

  The number of measurement cycles recording metabolic rate in an animal
  (a single number value 1 through 4).

- cycle_start:

  A numeric vector of measurement cycle start times (min). The length of
  this vector has to be equal to the number of cycles provided, e.g. if
  cycles = 2, then it may be cycle_start = c(0,10))

- cycle_end:

  A numeric vector of measurement cycle end times (min). The same format
  as cycle_start

- mmr_Ch1:

  A number indicated in which cycle the animal in oxygen meter Channel 1
  is expected to reach its maximum metabolic rate. Must specify to
  analyze metabolic rate, or enter 0, which will drop the channel from
  the analysis.

- mmr_Ch2:

  See mmr_Ch1 argument for description. This applies to Channel 2

- mmr_Ch3:

  See mmr_Ch1 argument for description. This applies to Channel 3

- mmr_Ch4:

  See mmr_Ch1 argument for description. This applies to Channel 4

- clean_Ch1:

  Specific to MMR cycle for each animal (or cycle). Use to set a time
  window (min) that specifies a quality section on metabolic rate
  measurement for sliding window analysis. Default will use full length
  measurement cycle as set by cycle_start and cycle_end arguments.

- clean_Ch2:

  See mmr_Ch1 argument for description. This applies to Channel 2

- clean_Ch3:

  See mmr_Ch1 argument for description. This applies to Channel 3

- clean_Ch4:

  See mmr_Ch1 argument for description. This applies to Channel 4

- N_Ch:

  The number of oxygen and temperature channels. Options include 2, 4,
  and 8. This argument can be ignored if 2-Channel oxygen meter was
  used.

- local_path:

  Logical. If TRUE (default) all returned files will be saved in the
  local working directory.

- date_format:

  The date format used in the original files. Must specify one of the
  following: c("m/d/y","d/m/y","y-m-d")

- inventory_data:

  The name of an inventory data dataframe that provides details of
  "cleaning" requirements

## Value

The output from [`print`](https://rdrr.io/r/base/print.html)
