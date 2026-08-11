# Estimates regressions for automatically timed measurement cycles

Analyzes O~2~ and temperature data recorded continuously using repeat
open(flush):closed(measure) cycles to obtain regressions representative
of metabolic rates.

## Usage

``` r
SMR(
  data,
  cycle_start,
  cycle_end,
  chop_start,
  chop_end,
  first_cycle = "flush",
  N_Ch = 4,
  date_format = c("%m/%d/%Y %H:%M:%S", "GMT"),
  inventory_data = NULL,
  local_path = TRUE,
  plot_temp = FALSE,
  background_data = FALSE,
  sda_data = FALSE,
  name_extensionID = "",
  drop_ch = 0
)
```

## Arguments

- data:

  The name of the raw input file .csv, character.

- cycle_start:

  The stat time (min) of auto cycle. This indicates the start of the
  measurement (relative to min 0).

- cycle_end:

  The time (min) of auto cycle; end of the measurement (relative to min
  0)

- chop_start:

  The time (min) to be chopped off at the start of each measurement
  cycle

- chop_end:

  The time (min) to be chopped off at the end of each measurement cycle

- first_cycle:

  whether the first cycle is a flush or metabolic rate measurement.
  provide either: "flush" or "measure"

- N_Ch:

  The number of channels of the oxygen meter. It must be either 4 or 8.
  8-Channel oxygen meter has 4 temperature and 4 oxygen sensors. The
  microplate reader files have 4 at this step.

- date_format:

  Common format of from the CSV conversion.

- inventory_data:

  The inventory file (.csv file), default = NULL

- local_path:

  Logical. If TRUE (default) all returned files will be saved in the
  local working directory.

- plot_temp:

  Logical argument. Indicates whether or not temperature trends for each
  cycle will be plotted and saved TRUE

- background_data:

  logical. If this datafile is a background (background respiration
  file), indicate TRUE

- sda_data:

  logical. If this datafile belong to SDA analysis, is a SDA data file,
  indicate TRUE

- name_extensionID:

  add identifier to the original name of the file

- drop_ch:

  exclude channel from analysis; default = 0 (none)

## Value

The output from [`print`](https://rdrr.io/r/base/print.html)
