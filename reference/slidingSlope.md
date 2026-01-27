# steepest slope analysis only

Estimates slopes at every time iteration (often second), and pulls out
the steepest 30 s, 60 s, 90 s, and 180 s slopes, as well as the
regression slope for the entire recorded duration.

## Usage

``` r
slidingSlope(data, Ch, local_path, N_Ch)
```

## Arguments

- data:

  Dataframe with each measurement present (O2 observed each second )

- Ch:

  Channel that needs to be analyzed (use format: c(1, 2), c(3), etc. up
  to 4 channels)

- local_path:

  Logical. If TRUE (no default) all returned files will be saved in the
  local working directory. Can also provide a path if this function is
  run independently

- N_Ch:

  Number of channels for the oxygen meter

## Value

The output from [`print`](https://rdrr.io/r/base/print.html)
