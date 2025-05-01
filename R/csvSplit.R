
#' Function to split to csv files for independent analyses.
#'
#' @param data The data file that needs to be split (this is after the raw file .txt/.xlsx/.csv conversion to data-analysis csv file)
#' @param split_data_name Some unique names string for the new files (two files split at the indicated time)
#' @param time_split Time at which the split should happen (in minutes)
#' @param date_format Common format of from the CSV conversion.
#' @param split Logical, if TRUE, the two new files will be saved; if FALSE the plot for visualization of the split will be saved, but the data will not.
#' @param N_Ch The number of channels in the csv file. It must be either 4 or 8. 8-Channel oxygen meter has 4 temperature and 4 oxygen sensors. The microplate reader files have 4 at this step.
#'
#' @return The output from \code{\link{print}}
#' @export
#'
csvSplit <- function (data,
                      split_data_name,
                      time_split,
                      date_format = c("%m-%d-%Y %H:%M:%S", "GMT"),
                      split=FALSE,
                      N_Ch=4 ) {

  if(!length(as.vector(date_format))==2){
    stop_function<-TRUE
    if(stop_function) {
      stop("Argument 'date_format' is not properly stated. \n It must be a vector of two stating: i) the date time format and ii) timezone. \n Default is: c(\"%m/%d/%Y %H:%M:%S\", \"GMT\"). Argument is passed to strptime() ")
    }
  }

  message("Output files saved locally")

  # names for the export files
  plotname1<-paste(split_data_name, "_split_respo_cycles.png", sep="")
  plotname2<-paste(split_data_name, "_split_respo_confirm.png", sep="")

  # read in file adn
  data_full<-read.csv(data)

  # print(nrow(new_csv))
  data_full<-data_full[!(is.na(data_full$Ch1_O2) &
          is.na(data_full$Ch2_O2) &
          is.na(data_full$Ch3_O2) &
          is.na(data_full$Ch4_O2)),] # all channel o2 is NA

  data_full$date<-as.character(data_full$date)
  data_full$time<-as.character(data_full$time)

	if(grepl(pattern = "M", as.character(data_full$time[1]))){
    data_full$time<-format(strptime(data_full$time, format = '%I:%M:%S %p'), format='%H:%M:%S')
	}

  # format time
	DateTime<-strptime(paste(data_full$date, data_full$time), format = date_format[1], tz = date_format[2])

  if(is.na(DateTime[1])){
    message(paste("DateTime is NA, likely wrong format provided: is it", date_format[1], "?"))
  }

  data_full$DateTime<-DateTime
  data_full$time_sec2[2:nrow(data_full)]<-diff(data_full$time_sec)
  data_full$hr<-round(data_full$time_sec/3600,2)
  data_full$hr2<-round(data_full$time_sec/3600,0)
  data_full$time_min<-round(data_full$time_sec/60,2)


	if(c(as.character(data_full$Ch1_O2[1])=="--- " )||
	   c(as.character(data_full$Ch1_O2[1])=="---" )||
	   c(is.na(data_full$Ch1_O2[1]) )){
		data_full$Ch1_O2<-0
	}
	if(c(as.character(data_full$Ch2_O2[1])=="--- ")||
	   c(as.character(data_full$Ch2_O2[1])=="---" )||
	   c(is.na(data_full$Ch2_O2[1]) )){
		data_full$Ch2_O2<-0
	}
	if(c(as.character(data_full$Ch3_O2[1])=="--- " )||
	   c(as.character(data_full$Ch3_O2[1])=="---" )||
	   c(is.na(data_full$Ch3_O2[1]) )){
		data_full$Ch3_O2<-0
	}
	if(c(as.character(data_full$Ch4_O2[1])=="--- " )||
	   c(as.character(data_full$Ch4_O2[1])=="---" )||
	   c(is.na(data_full$Ch4_O2[1]) )){
		data_full$Ch4_O2<-0
	}

  # times in minutes
  time_split0 <- time_split-data_full$time_min[2]-data_full$time_min[1]

  # rows with correct times
  # rows to split at; last row for file 1

  t0_row <- which(abs(data_full$time_min - time_split0) == min(abs(data_full$time_min - time_split0), na.rm = TRUE))[1]

  # the next row for file 1
  t1_row <- which(abs(data_full$time_min - time_split) == min(abs(data_full$time_min - time_split), na.rm = TRUE))[1]

  png(plotname1, width = 10, height = 20, units="in", res=200)
    par(mfrow=c(4,1))
    plot(data_full$time_min, data_full$Ch1_O2, main="Channel 1",
         cex = 0.1, xlab = "Time (min)", ylab = "Oxygen (mgO2/L)")
    abline(v=data_full$time_min[t0_row], col="blue", lwd=1)

    plot(data_full$time_min, data_full$Ch2_O2, main="Channel 2",
         cex = 0.1, xlab = "Time (min)", ylab = "Oxygen (mgO2/L)")
    abline(v=data_full$time_min[t0_row], col="blue", lwd=1)

    plot(data_full$time_min, data_full$Ch3_O2, main="Channel 3",
         cex = 0.1, xlab = "Time (min)", ylab = "Oxygen (mgO2/L)")
    abline(v=data_full$time_min[t0_row], col="blue", lwd=1)

    plot(data_full$time_min, data_full$Ch4_O2, main="Channel 4",
         cex = 0.1, xlab = "Time (min)", ylab = "Oxygen (mgO2/L)")
    abline(v=data_full$time_min[t0_row], col="blue", lwd=1)
  dev.off()

  # full file splitting part
  if(split==TRUE){

  	file1name<-paste(split_data_name, "_split_file1.csv", sep="")
    file2name<-paste(split_data_name, "_split_file2.csv", sep="")

  	data1<-data_full[1:t0_row,]
  	data2<-data_full[t1_row:nrow(data_full),]

  	data2$time_sec<-abs(data2$time_sec[1] - data2$time_sec)

    data1$time_min<-data1$time_sec/60
  	data2$time_min<-data2$time_sec/60

    # get rid of extra columns
    if(N_Ch != 4){
      data1<-data1[,1:11]
      data2<-data2[,1:11]
    }else{
      data1<-data1[,1:8]
      data2<-data2[,1:8]
    }
  	  write.csv(file=file1name, data1, row.names=FALSE)
  	  write.csv(file=file2name, data2, row.names=FALSE)
  }else{
    message("Data outputs not saved. To save add/change argument: split = TRUE")
  }

}
