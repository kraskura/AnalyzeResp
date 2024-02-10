#' Title
#'
#' @param txt_file The name of the original “.txt” file from FireSting
#' @param nrowSkip The number of rows to skip if the default does not work; the argument 'type_file' determines the default N rows to skip "Firesting_pre2023" = 19, "Firesting_2023" = 70, "Witrox" = 41; these rows in raw txt file often contain calibration information, the IDs of the probes and other user defined settings
#' @param type_file Indicates the type of software that was used to record raw data, options: "Firesting_pre2023", "Firesting_2023", "Witrox"
#' @param N_Ch The number of FireSting channels. Options include 2, 4, 8. If a 2-channel FireSting was used, this argument could be ignored, or enter 4
#' @param local_path Logical. If TRUE (default) all returned files will be saved in the local working directory.
#' @param exclude_first_measurement_s The number measurement point to be excluded from the beginning of the file (in addition to the nrowSkip argument)
#' @param convert_units Logical (FALSE = default). If true, the function is passed to rMR function DO.unit.convert to convert O2 content units.
#' @param units_from default NULL, options:"mg/L", "PP", "pct". passed down to rM::DO.unit.convert arg. DO.units.in
#' @param units_to = default NULL, options:"mg/L", "PP", "pct". passed down to rM::DO.unit.convert arg. DO.units.out
#' @param channels = c(1,2,3,4), indicate which channels to apply the conversion to
#' @param salinity = 0, passed down to rM::DO.unit.convert arg. salinity (must be in "pp.thou")
#' @param atm_pressure = 1, passed down to rM::DO.unit.convert arg. bar.press (must be in "atm")
#' @param temperature user set consistent temperature for all channels (ºC)
#'
#' @return The output from \code{\link{print}}
#' @export
#'
#' @importFrom utils read.delim
#' @importFrom rMR DO.unit.convert
#'
textFileConvert<-function(txt_file,
                          type_file,
                          nrowSkip = NULL,
                          N_Ch = 4,
                          local_path = TRUE,
                          exclude_first_measurement_s = 0,
                          convert_units = FALSE,
                          units_from = NULL,
                          units_to = NULL,
                          channels = c(1,2,3,4),
                          salinity = 0,
                          atm_pressure = 1,
                          temperature = NULL){

  # if(!is.numeric(nrowSkip)){
  #   stop("Must provide how many rows to skip from raw datafile;
  #        e.g. 4-Ch and 2-Ch fireSting commonly need 19, 8-Ch Firesting needs 26")
  # }

  if(type_file == "Firesting_pre2023"){

    if(N_Ch == 4 | N_Ch==2){
    	new_csv<-as.data.frame(matrix(nrow=0, ncol=8))
    	colnames(new_csv)<-c("date", "time", "time_sec", "Ch1_O2", "Ch1_temp", "Ch2_O2", "Ch3_O2", "Ch4_O2")
      if(is.null(nrowSkip)){
        nrowSkip <- 19
      }
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

    	if(is.null(temperature)){
    	   new_csv$Ch1_temp<-d[,15]# temp Ch1 - but same for all
    	}else{
         new_csv$Ch1_temp<-temperature
    	}

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

    	if(is.null(temperature)){
      	new_csv$Ch1_temp<-d[,9]# unique ch temps
      	new_csv$Ch2_temp<-d[,10]# unique ch temps
      	new_csv$Ch3_temp<-d[,11]# unique ch temps
      	new_csv$Ch4_temp<-d[,12]# unique ch temps
    	}else{
      	new_csv$Ch1_temp<-temperature
      	new_csv$Ch2_temp<-temperature
      	new_csv$Ch3_temp<-temperature
      	new_csv$Ch4_temp<-temperature
    	}

    	new_csv$Ch2_O2<-d[,6]
    	new_csv$Ch3_O2<-d[,7]
    	new_csv$Ch4_O2<-d[,8]
    }

  }
  else if (type_file == "Firesting_2023"){
      new_csv<-as.data.frame(matrix(nrow=0, ncol=8))
      if(is.null(nrowSkip)){
        nrowSkip <-70
      }
    	d<-read.delim(txt_file, skip = nrowSkip + exclude_first_measurement_s) # 70
      colnames(d)<-d[1,]
      d<-d[-1,]

    	names<-colnames(d)
      O2_ch1_name<-which(c(grepl("Oxygen", x = names, ignore.case = T) & grepl("Ch.1", x = names, ignore.case = T)))
      O2_ch2_name<-which(c(grepl("Oxygen", x = names, ignore.case = T) & grepl("Ch.2", x = names, ignore.case = T)))
      O2_ch3_name<-which(c(grepl("Oxygen", x = names, ignore.case = T) & grepl("Ch.3", x = names, ignore.case = T)))
      O2_ch4_name<-which(c(grepl("Oxygen", x = names, ignore.case = T) & grepl("Ch.4", x = names, ignore.case = T)))
      temp_ch1_name<-which(c(grepl("temp", x = names, ignore.case = T) & grepl("Ch.1", x = names, ignore.case = T)))
      temp_ch2_name<-which(c(grepl("temp", x = names, ignore.case = T) & grepl("Ch.2", x = names, ignore.case = T)))
      temp_ch3_name<-which(c(grepl("temp", x = names, ignore.case = T) & grepl("Ch.3", x = names, ignore.case = T)))
      temp_ch4_name<-which(c(grepl("temp", x = names, ignore.case = T) & grepl("Ch.4", x = names, ignore.case = T)))    	# d<-d[,1:15]

      new_csv<-d[, c(1:3)]

    if(N_Ch == 4 | N_Ch==2){

    	new_csv$Ch1_O2<-d[,O2_ch1_name]
    	if(is.null(temperature)){
    	  new_csv$Ch1_temp<-d[,temp_ch1_name]# temp Ch1 - but same for all
    	}else{
    	  new_csv$Ch1_temp<-temperature
    	}

    	new_csv$Ch2_O2<-d[,O2_ch2_name]
    	new_csv$Ch3_O2<-d[,O2_ch3_name]
    	new_csv$Ch4_O2<-d[,O2_ch4_name]

    	colnames(new_csv)<-c("date", "time", "time_sec", "Ch1_O2", "Ch1_temp", "Ch2_O2", "Ch3_O2", "Ch4_O2")
    }

    if(N_Ch==8){
      if(is.null(temperature)){
      	new_csv$Ch2_temp<-d[,temp_ch2_name]
        new_csv$Ch3_temp<-d[,temp_ch3_name]
        new_csv$Ch4_temp<-d[,temp_ch4_name]
      }else{
      	new_csv$Ch2_temp<-temperature
        new_csv$Ch3_temp<-temperature
        new_csv$Ch4_temp<-temperature
      }

    	colnames(new_csv)<-c("date", "time", "time_sec", "Ch1_O2", "Ch1_temp", "Ch2_O2", "Ch3_O2", "Ch4_O2", "Ch2_temp", "Ch3_temp", "Ch4_temp")
    }

  }
  else if (type_file == "Witrox_2023"){
      new_csv<-as.data.frame(matrix(nrow=0, ncol=8))
      if(is.null(nrowSkip)){
        nrowSkip <-41
      }
      d<-read.table(txt_file, skip = nrowSkip + exclude_first_measurement_s,
                    sep = "\t", skipNul = TRUE, header = FALSE)
      colnames(d)<-d[1,]
      d<-d[-1,]

    	names<-colnames(d)
      O2_ch1_name<-which(c(grepl("Oxygen", x = names, ignore.case = T) & grepl("Ch 1", x = names, ignore.case = T)))
      O2_ch2_name<-which(c(grepl("Oxygen", x = names, ignore.case = T) & grepl("Ch 2", x = names, ignore.case = T)))
      O2_ch3_name<-which(c(grepl("Oxygen", x = names, ignore.case = T) & grepl("Ch 3", x = names, ignore.case = T)))
      O2_ch4_name<-which(c(grepl("Oxygen", x = names, ignore.case = T) & grepl("Ch 4", x = names, ignore.case = T)))
      temp_ch1_name<-which(c(grepl("temp", x = names, ignore.case = T) & grepl("Ch 1", x = names, ignore.case = T)))
      temp_ch2_name<-which(c(grepl("temp", x = names, ignore.case = T) & grepl("Ch 2", x = names, ignore.case = T)))
      temp_ch3_name<-which(c(grepl("temp", x = names, ignore.case = T) & grepl("Ch 3", x = names, ignore.case = T)))
      temp_ch4_name<-which(c(grepl("temp", x = names, ignore.case = T) & grepl("Ch 4", x = names, ignore.case = T)))    	# d<-d[,1:15]

      new_csv<-d[, c(1:2)]

      new_csv$time_sec<-
        sapply(strsplit(d$`Relative time [HH:MM:SS]`,":"),
        function(x) {
          x <- as.numeric(x)
          x[3]+x[2]*60+x[1]*60*60
          }
      )

    	new_csv$Ch1_O2<-d[,O2_ch1_name]
    	new_csv$Ch1_temp<-d[,temp_ch1_name]# temp Ch1 - but same for all

    	new_csv$Ch2_O2<-d[,O2_ch2_name]
    	new_csv$Ch3_O2<-d[,O2_ch3_name]
    	new_csv$Ch4_O2<-d[,O2_ch4_name]

    	colnames(new_csv)<-c("date", "time", "time_sec", "Ch1_O2", "Ch1_temp", "Ch2_O2", "Ch3_O2", "Ch4_O2")
    }

  new_csv[c(3:ncol(new_csv))] <- sapply(new_csv[3:ncol(new_csv)],as.numeric)
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
      new_csv$Ch1_O2 <- DO.unit.convert(new_csv$Ch1_O2, DO.units.in = units_from,
                                        DO.units.out = units_to, bar.units.in ="atm",
                    bar.press = atm_pressure, temp.C = new_csv$Ch1_temp, bar.units.out = "atm",
                    salinity = salinity, salinity.units = "pp.thou")

    }

    if(any(channels == 2)){
      # new_csv$Ch2_O2 <- conv_o2 (o2 = new_csv$Ch2_O2, from = "percent_a.s.", to = "mg_per_l", temp = temp_ch2, sal = salinity, atm_pres = atm_pressure)
      if (N_Ch!= 8){
          new_csv$Ch2_O2 <- DO.unit.convert(new_csv$Ch2_O2, DO.units.in = units_from, DO.units.out = units_to, bar.units.in ="atm",
                      bar.press = atm_pressure, temp.C = new_csv$Ch1_temp, bar.units.out = "atm",
                      salinity = salinity, salinity.units = "pp.thou")

      }else{
          new_csv$Ch2_O2 <- DO.unit.convert(new_csv$Ch2_O2, DO.units.in = units_from, DO.units.out = units_to, bar.units.in ="atm",
                      bar.press = atm_pressure, temp.C = new_csv$Ch2_temp, bar.units.out = "atm",
                      salinity = salinity, salinity.units = "pp.thou")
      }

    }

    if(any(channels == 3)){
      # new_csv$Ch3_O2 <- conv_o2 (o2 = new_csv$Ch3_O2, from = "percent_a.s.", to = "mg_per_l", temp = temp_ch3, sal = salinity, atm_pres = atm_pressure)
      if (N_Ch!= 8){
        new_csv$Ch3_O2 <- DO.unit.convert(new_csv$Ch3_O2, DO.units.in = units_from, DO.units.out = units_to, bar.units.in ="atm",
                      bar.press = atm_pressure, temp.C = new_csv$Ch1_temp, bar.units.out = "atm",
                      salinity = salinity, salinity.units = "pp.thou")
      }else{
        new_csv$Ch3_O2 <- DO.unit.convert(new_csv$Ch3_O2, DO.units.in = units_from, DO.units.out = units_to, bar.units.in ="atm",
                                                    bar.press = atm_pressure, temp.C = new_csv$Ch3_temp, bar.units.out = "atm",
                                                    salinity = salinity, salinity.units = "pp.thou")
      }
    }

    if(any(channels == 4)){
      # new_csv$Ch4_O2 <- conv_o2 (o2 = new_csv$Ch4_O2, from = "percent_a.s.", to = "mg_per_l", temp = temp_ch4, sal = salinity, atm_pres = atm_pressure)
      if (N_Ch!= 8){
        new_csv$Ch4_O2 <- DO.unit.convert(new_csv$Ch4_O2, DO.units.in = units_from, DO.units.out = units_to, bar.units.in ="atm",
                      bar.press = atm_pressure, temp.C = new_csv$Ch1_temp, bar.units.out = "atm",
                      salinity = salinity, salinity.units = "pp.thou")
      }else{
        new_csv$Ch4_O2 <- DO.unit.convert(new_csv$Ch4_O2, DO.units.in = units_from, DO.units.out = units_to, bar.units.in ="atm",
                                                    bar.press = atm_pressure, temp.C = new_csv$Ch4_temp, bar.units.out = "atm",
                                                    salinity = salinity, salinity.units = "pp.thou")
      }
    }
  }

	# save in the current directory (default)
  if(local_path){
    write.csv(file=paste(gsub('.{4}$', '', txt_file), "_converted.csv", sep=''), new_csv, row.names=FALSE)
    message("Return csv files saved in local working directory")
  } else if(local_path == FALSE & dir.exists("csv_files")){
    write.csv(file=paste("./csv_files/", gsub('.{4}$', '', txt_file), "_converted.csv", sep=''), new_csv, row.names=FALSE)
	  message("Return csv files saved in \"./csv_files\" local directory")
	}else{
    stop("Cannot save new file, need either csv_files folder, or local_path = TRUE argument")
	}

}
