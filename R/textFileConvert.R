#' Title
#'
#' @param txt_file The name of the original “.txt” file from FireSting
#' @param nrowSkip The number of rows to skip; these rows in raw txt file often contain calibration information, the IDs of the probes and other user defined settings
#' @param N_Ch The number of FireSting channels. Options include 2, 4, 8. If a 2-channel FireSting was used, this argument could be ignored, or enter 4
#' @param local_path Logical. If TRUE (default) all returned files will be saved in the local working directory.
#' @param exclude_first_measurement_s The number measurement point to be excluded from the beginning of the file (in addition to the nrowSkip argument)
#'
#' @return The output from \code{\link{print}}
#' @export
#'
#' @importFrom utils read.delim
#'
textFileConvert<-function(txt_file,
                          nrowSkip,
                          N_Ch = 4,
                          local_path = TRUE,
                          exclude_first_measurement_s = 0,
                          convert_units = FALSE,
                          units_from = NULL,
                          units_to = NULL,
                          channels = c(1,2,3,4),
                          salinity = 0,
                          atm_pressure = 1){

  if(!is.numeric(nrowSkip)){
    stop("Must provide how many rows to skip from raw datafile; e.g. 4-Ch and 2-Ch fireSting commonly need 19, 8-Ch Firesting needs 26")
  }

  if(N_Ch == 4 | N_Ch==2){
  	new_csv<-as.data.frame(matrix(nrow=0, ncol=8))
  	colnames(new_csv)<-c("date", "time", "time_sec", "Ch1_O2", "Ch1_temp", "Ch2_O2", "Ch3_O2", "Ch4_O2")

  	d<-read.delim(txt_file, skip = nrowSkip + exclude_first_measurement_s)

  	d<-d[,1:15]

  	nr<-nrow(d)
  	nc<-ncol(d)

  	new_csv[nr,]<-NA
  	new_csv$date<-d[,1]
  	new_csv$time<-d[,2]
  	new_csv$time_sec<-d[,3]
  	# oxyegen below
  	new_csv$Ch1_O2<-d[,5]

  	new_csv$Ch1_temp<-d[,15]# temp Ch1 - but same for all

  	new_csv$Ch2_O2<-d[,6]
  	new_csv$Ch3_O2<-d[,7]
  	new_csv$Ch4_O2<-d[,8]
  } # end od N_Ch == 2

  if(N_Ch==8){
    new_csv<-as.data.frame(matrix(nrow=0, ncol=11))
  	colnames(new_csv)<-c("date", "time", "time_sec", "Ch1_O2", "Ch1_temp", "Ch2_O2", "Ch3_O2", "Ch4_O2", "Ch2_temp", "Ch3_temp", "Ch4_temp")

  	d<-read.delim(txt_file, skip= nrowSkip + exclude_first_measurement_s)
  	# exclude rows with wrong data, no data, no time formats for the two columns
  	d<-d[c(which(grepl( "/", d[,1]) & grepl( ":", d[,2]))),]
  	d<-d[,1:12]
  	nr<-nrow(d)
  	nc<-ncol(d)

  	new_csv[nr,]<-NA
  	new_csv$date<-d[,1]
  	new_csv$time<-d[,2]
  	new_csv$time_sec<-d[,3]
  	# oxyegen below
  	new_csv$Ch1_O2<-d[,5]

  	new_csv$Ch1_temp<-d[,9]# unique ch temps
  	new_csv$Ch2_temp<-d[,10]# unique ch temps
  	new_csv$Ch3_temp<-d[,11]# unique ch temps
  	new_csv$Ch4_temp<-d[,12]# unique ch temps

  	new_csv$Ch2_O2<-d[,6]
  	new_csv$Ch3_O2<-d[,7]
  	new_csv$Ch4_O2<-d[,8]
  }

  if(N_Ch==4){
    temp_ch1 <- new_csv$Ch1_temp
    temp_ch2 <- new_csv$Ch1_temp
    temp_ch3 <- new_csv$Ch1_temp
    temp_ch4 <- new_csv$Ch1_temp
  }
  if(N_Ch==2){
    temp_ch1 <- new_csv$Ch1_temp
    temp_ch2 <- new_csv$Ch1_temp
  }
  if(N_Ch==8){
    temp_ch1 <- new_csv$Ch1_temp
    temp_ch2 <- new_csv$Ch2_temp
    temp_ch3 <- new_csv$Ch3_temp
    temp_ch4 <- new_csv$Ch4_temp
  }

  if(convert_units){
    if(is.null(units_from) | is.null(units_to)){
      stop_function<-TRUE
      if(stop_function){
        stop("If converting units, must provide 'units_from' and 'units_to'")
      }
    }
    message(paste("Unit conversion parameters:","\n",
            " units_from: ", units_from, "\n",
            " units_to: ", units_to, "\n",
            " salinity: ", salinity, "\n",
            " atm_pressure: ", atm_pressure, sep =""))

    if(any(channels == 1)){
      # new_csv$Ch1_O2 <- conv_o2 (o2 = new_csv$Ch1_O2, from = "percent_a.s.", to = "mg_per_l", temp = temp_ch1, sal = salinity, atm_pres = atm_pressure)

      new_csv$Ch1_O2 <- rMR::DO.unit.convert(new_csv$Ch1_O2, DO.units.in = units_from, DO.units.out = units_to, bar.units.in ="atm",
                    bar.press=1, temp.C = temp_ch1, bar.units.out = "atm",
                    salinity = salinity, salinity.units = "pp.thou")

    }

    if(any(channels == 2)){
      # new_csv$Ch2_O2 <- conv_o2 (o2 = new_csv$Ch2_O2, from = "percent_a.s.", to = "mg_per_l", temp = temp_ch2, sal = salinity, atm_pres = atm_pressure)
      if (N_Ch!= 8){
          new_csv$Ch2_O2 <-  rMR::DO.unit.convert(new_csv$Ch2_O2, DO.units.in = units_from, DO.units.out = units_to, bar.units.in ="atm",
                      bar.press=1, temp.C = temp_ch1, bar.units.out = "atm",
                      salinity = salinity, salinity.units = "pp.thou")

      }else{
          new_csv$Ch2_O2 <- rMR::DO.unit.convert(new_csv$Ch2_O2, DO.units.in = units_from, DO.units.out = units_to, bar.units.in ="atm",
                      bar.press=1, temp.C = temp_ch2, bar.units.out = "atm",
                      salinity = salinity, salinity.units = "pp.thou")
      }

    }

    if(any(channels == 3)){
      # new_csv$Ch3_O2 <- conv_o2 (o2 = new_csv$Ch3_O2, from = "percent_a.s.", to = "mg_per_l", temp = temp_ch3, sal = salinity, atm_pres = atm_pressure)
      if (N_Ch!= 8){
        new_csv$Ch3_O2 <- rMR::DO.unit.convert(new_csv$Ch3_O2, DO.units.in = units_from, DO.units.out = units_to, bar.units.in ="atm",
                      bar.press=1, temp.C = temp_ch1, bar.units.out = "atm",
                      salinity = salinity, salinity.units = "pp.thou")
      }else{
        new_csv$Ch3_O2 <- rMR::DO.unit.convert(new_csv$Ch3_O2, DO.units.in = units_from, DO.units.out = units_to, bar.units.in ="atm",
                                                    bar.press=1, temp.C = temp_ch3, bar.units.out = "atm",
                                                    salinity = salinity, salinity.units = "pp.thou")
      }
    }

    if(any(channels == 4)){
      # print(new_csv)
      # new_csv$Ch4_O2 <- conv_o2 (o2 = new_csv$Ch4_O2, from = "percent_a.s.", to = "mg_per_l", temp = temp_ch4, sal = salinity, atm_pres = atm_pressure)
      if (N_Ch!= 8){
        new_csv$Ch4_O2 <- rMR::DO.unit.convert(new_csv$Ch4_O2, DO.units.in = units_from, DO.units.out = units_to, bar.units.in ="atm",
                      bar.press=1, temp.C = temp_ch1, bar.units.out = "atm",
                      salinity = salinity, salinity.units = "pp.thou")
      }else{
        new_csv$Ch4_O2 <- rMR::DO.unit.convert(new_csv$Ch4_O2, DO.units.in = units_from, DO.units.out = units_to, bar.units.in ="atm",
                                                    bar.press=1, temp.C = temp_ch4, bar.units.out = "atm",
                                                    salinity = salinity, salinity.units = "pp.thou")
      }
    }
  }

	# save in the current directory (default)
  if(local_path | !dir.exists("csv_files")){
    write.csv(file=paste(gsub('.{4}$', '', txt_file), ".csv", sep=''), new_csv, row.names=FALSE)
    if(!dir.exists("csv_files")){
      message("Return csv files saved in local working directory")
    }
  }

	if(dir.exists("csv_files")){
    write.csv(file=paste("./csv_files/", gsub('.{4}$', '', txt_file), ".csv", sep=''), new_csv, row.names=FALSE)
	 message("Return csv files saved in \"./csv_files\" local directory")
	}

}
