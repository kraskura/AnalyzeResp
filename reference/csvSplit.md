# Function to split to csv files for independent analyses.

Function to split to csv files for independent analyses.

## Usage

``` r
csvSplit(
  data,
  split_data_name,
  time_split,
  date_format = c("%m-%d-%Y %H:%M:%S", "GMT"),
  split = FALSE,
  N_Ch = 4
)
```

## Arguments

- data:

  The data file that needs to be split (this is after the raw file
  .txt/.xlsx/.csv conversion to data-analysis csv file)

- split_data_name:

  Some unique names string for the new files (two files split at the
  indicated time)

- time_split:

  Time at which the split should happen (in minutes)

- date_format:

  Common format of from the CSV conversion.

- split:

  Logical, if TRUE, the two new files will be saved; if FALSE the plot
  for visualization of the split will be saved, but the data will not.

- N_Ch:

  The number of channels in the csv file. It must be either 4 or 8.
  8-Channel oxygen meter has 4 temperature and 4 oxygen sensors. The
  microplate reader files have 4 at this step.

## Value

The output from [`print`](https://rdrr.io/r/base/print.html)
