# filters and converts raw data files

Uses to convert and format (strip from unused columns and rows) the raw
'.txt' data file to produce a .csv file that is required for \`MMR\` and
\`SMR\` functions.

## Usage

``` r
textFileConvert(
  txt_file,
  type_file,
  nrowSkip = NULL,
  N_Ch = 4,
  local_path = TRUE,
  exclude_rows = NULL,
  exclude_first_measurement_s = 0,
  convert_units = FALSE,
  units_from = NULL,
  units_to = NULL,
  channels = c(1, 2, 3, 4),
  salinity = 0,
  atm_pressure = 1,
  temperature = NULL,
  temperature_Ch = 1,
  device = "",
  file_extension_id = ""
)
```

## Arguments

- txt_file:

  The name of the original “.txt” file from oxygen meter

- type_file:

  Indicates the type of software that was used to record raw data,
  options: "Firesting_pre2023", "Firesting_2023", "Witrox_pre2023",
  "Witrox_2023", "PreSense_microplate"

- nrowSkip:

  The number of rows to skip if the default does not work; the argument
  'type_file' determines the default N rows to skip "Firesting_pre2023"
  = 19, "Firesting_2023" = 70, "Witrox_2023" = 41; "PreSense_microplate"
  = 10, these rows in raw txt file often contain calibration
  information, the IDs of the probes and other user defined settings

- N_Ch:

  The number of oxygen meter channels that described the physical
  device, not how many channels were plugged in; this is an argument
  that describes hardware. Options include 2, 4, 8, 24. (microplate). If
  a 2-channel oxygen meter was used, this argument could be ignored, or
  enter 4

- local_path:

  Logical. If TRUE (default) all returned files will be saved in the
  local working directory.

- exclude_rows:

  Rows to be excluded; used when there are unwanted NAs, or no sensor
  data. etc.

- exclude_first_measurement_s:

  DEPRECATED Nov 2024: use 'exclude_measurement_s' (see documentation),
  The number measurement point to be excluded from the beginning of the
  file (in addition to the nrowSkip argument)

- convert_units:

  Logical (FALSE = default). If true, the function is passed to rMR
  function DO.unit.convert to convert O2 content units.

- units_from:

  default NULL, options:"mg/L", "PP", "pct", "mmol/L". passed down to
  rM::DO.unit.convert arg. DO.units.in. Note that "mmol/L" can be only
  used to go to "mg/L".

- units_to:

  = default NULL, options:"mg/L", "PP", "pct". passed down to
  rM::DO.unit.convert arg. DO.units.out

- channels:

  = c(1,2,3,4), indicate which channels the unit conversion will be
  applied to

- salinity:

  = 0, passed down to rMR::DO.unit.convert arg. salinity (must be in
  "pp.thou")

- atm_pressure:

  = 1, passed down to rMR::DO.unit.convert arg. bar.press (must be in
  "atm")

- temperature:

  user set consistent temperature for all channels (ºC)

- temperature_Ch:

  Numerical, select the channel to use for temperature recording. This
  is used only for 4-channel Firesting, when probe/channel 1 was not
  plugged in and therefore has all values as "NA". default is 'NULL'0'
  which then chooses data from the first available channel

- device:

  only for for "Firesting_2023" output files. Results from multiple
  devices can be recorded on one file. These are differentiated by
  uppper case letters "A", "B", etc. use this argument to specify which
  device is used, default is "A".

- file_extension_id:

  custom file extension for the saved csv file, default is 'blank'

## Value

The output from [`print`](https://rdrr.io/r/base/print.html)
