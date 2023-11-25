
#' Title
#'
#' @param data The name of the raw input file .csv, character.
#' @param cycle_start The stat time (min) of auto cycle. This indicates the start of the measurement (relative to min 0).
#' @param cycle_end The time (min) of auto cycle; end of the measurement (relative to min 0)
#' @param first_cycle whether the first cycle is a flush or metabolic rate measurement. provide either: "flush" or "measurement"
#' @param chop_start The time (min) to be chopped off at the start of each measurement cycle
#' @param chop_end The time (min) to be chopped off at the end of each measurement cycle
#' @param N_Ch The number of channels of the oxygen meter. It must be either 4 or 8. 8-Channel Firesting has 4 temperature and 4 oxygen sensors
#' @param inventory_data The inventory file (.csv file), default = NULL
#' @param local_path Logical. If TRUE (default) all returned files will be saved in the local working directory.
#' @param date_format Common format of from the CSV conversion.
#' @param plot_temp Logical argument. Indicates whether or not temperature trends for each cycle will be plotted and saved	TRUE
#' @param background_data logical. If this datafile is a background (background respiration file), indicate TRUE
#' @param sda_data logical. If this datafile belong to SDA analysis, is a SDA data file, indicate TRUE
#'
#' @return The output from \code{\link{print}}
#' @export
#'
#' @importFrom stats lm
#' @importFrom stats coef
#' @import graphics
#' @import grDevices
#' @import scales
#' @import ggplot2
#' @import utils
#'
#'
SMR<-function(data,
              cycle_start,
              cycle_end,
              first_cycle,
              chop_start,
              chop_end,
              N_Ch = 4,
              date_format = c("%m/%d/%Y %H:%M:%S", "GMT"),
              inventory_data = NULL,
              local_path = TRUE,
              plot_temp = FALSE,
              background_data = FALSE,
              sda_data = FALSE){

  if(!length(as.vector(date_format))==2){
    stop_function<-TRUE
    if(stop_function) {
      stop("Argument 'date_format' is not properly stated. \n It must be a vector of two stating: i) the date time format and ii) timezone. \n Default is: c(\"%m/%d/%Y %H:%M:%S\", \"GMT\"). Argument is passed to strptime() ")
    }
  }

  #  binding global variables locally to the function.
  min_start<-m<-Ch<-NULL

  # **********************************************************
  # only relevant inside the RMR repeat. leave them inside here

  # this function is called inside the MMR_SMR function to analyze data of each channel
  # don't need to worry about the variable inside this function, the MMR_SMR provides them all
  Channel<-function(Ch,
                    temp,
                    seq_st,
                    seq_end,
                    plotname,
                    data1,
                    chop_start,
                    chop_end,
                    inv.data,
                    newdata,
                    plot_temp,
                    N_Ch2){ # Ch = colum of the channel not the actual channel

    if(plot_temp==TRUE){
      # temp plot name
      sub_directory_name<-gsub(pattern="plots_channel", replacement= "plots_channel_temperature", x= plotname)
      plotname_temperature<-gsub(pattern=".png", replacement= "_temp.png", x= sub_directory_name)
    }

  	if(!N_Ch2==8){
  		if (Ch==4){
  		  temp<-5
  		}
  		if (Ch==6){
  		  temp<-5
  		}
  		if (Ch==7){
  		}
  		if (Ch==8){
  		  temp<-5
  		}
  	}else{
     if (Ch==4){
    		  temp<-5
    		}
    		if (Ch==6){
    		  temp<-9
    		}
    		if (Ch==7){
    		  temp<-10
    		}
    		if (Ch==8){
    		  temp<-11
    		}
  	  }

    # no inventory data-----
  	if(nrow(inv.data)==0){


  		# run through each "measure cycle" and get slope, intercept, r2. All channel cycles plotted together in one file
  	  for (i in 1:length(seq_end)){

        # ************************************************************
    	  cycle_time<-seq_end[2]-seq_end[1]
  	    end_fullCycle<-seq_end[i]
        start_fullCycle<-end_fullCycle-cycle_time

        start_NoChop<-seq_st[i]
  			end_NoChop<-seq_end[i]

  			start<-seq_st[i]+chop_start #
  			end<-seq_end[i]-chop_end # chopping off last x min

  			d_fullCycle<-data1[c(which(data1$time_min>start_fullCycle & data1$time_min<end_fullCycle)),]
  			d_NoChop<-data1[c(which(data1$time_min>start_NoChop & data1$time_min<end_NoChop)),] # from teh entire files, get only the data section we need
        d<-data1[c(which(data1$time_min>start & data1$time_min<end)),] # from teh entire files, get only the data section we need

  			if(nrow(d) == 0){
          next_function <- TRUE
          if(Ch==4){ # print this only for the first channel - assuming this is system wide technical error
            message(paste("! Missing data: some repeat cycles cannot be analyzed. Skip cycles: min", round(start), " to", round(end)))
          }
          next
  			}

  			DateTime_start<-as.character(d$DateTime[1])
  			lm_coef<-coef(lm(d[,Ch]~d$time_min)) # get linear regression fit
  			r2<-round(summary(lm(d[,Ch]~d$time_min))$r.squared,3) # get r2

  			m<-round(lm_coef[2],5) # get slope
  			b<-round(lm_coef[1],2) # get intercept

  				# if this is the first cycle for this fish - start a new data file]
#
#   				if(newdata[1,1]=="new"){
#   					newdata<-matrix(ncol=15,nrow=0)
#   					colnames(newdata)<-c("time_frame","min_start", "r2" ,"b", "m" , "t_min", "t_max", "t_mean","O2_min", "O2_max", "O2_mean", "Ch", "DateTime_start", "type", "n_min")
#   				}

  				# if this is the first cycle for this fish - start a new plot

  				if (i == 1){
  					if(length(seq_st)<=100){
  						png(plotname, width=40, height=40, units="in",  res=200)
  						par(mfrow=c(10,10)) # This will fit 100 plots.
  					}
  					if(length(seq_st)>100 & length(seq_st)<=150){
  						png(plotname, width=40, height=60, units="in",  res=200)
  						par(mfrow=c(15,10)) # This will fit 150 plots
  					}
  					if(length(seq_st)>150 & length(seq_st)<=225){
  						png(plotname, width=50, height=50, units="in",  res=200)
  						par(mfrow=c(15,15)) # This will fit 150 plots
  					}
  					if(length(seq_st)>=225 & length(seq_st)<=375){
    			  	png(plotname, width=50, height=70, units="in",  res=200)
    				  par(mfrow=c(15,25)) # This will fit 375 plots
    		    }
    		    if(length(seq_st)>375 ){
    			   	png(plotname, width=65, height=80, units="in",  res=200)
    			  	par(mfrow=c(18,35))
    			  	# print("here")
    		    }
  				}

  			# # plotting each segment with slope, r2, and equations on it
  			plot(d_fullCycle[,Ch]~time_min, data=d_fullCycle, ylab=expression(paste(O2~(mg~L^{-1}))),xlab="Time (min)", col = "grey", cex = 0.5, pch = 1)
        points(d_NoChop[,Ch]~time_min, data=d_NoChop, col = "blue", cex = 0.5, pch = 1)
  			points(d[,Ch]~time_min, data=d, col = "black", cex = 0.5, pch = 1)
        abline(v=start_NoChop, col="grey", lty = "dashed")
        abline(lm(d[,Ch]~d$time_min), col="red",lwd = 2)
  			mtext(bquote(y == .(lm_coef[2])*x + .(lm_coef[1])), adj=1, padj=0, cex=0.8, line=0) # display equation
  			mtext(bquote(italic(R)^2 == .(format(r2, digits = 3))),adj=1, padj=0, cex=0.8, line=1)
  			# the bwlow code for the frames


  				if(r2<0.9){
  					box(lty="dashed", col="orange", lwd=2)
  				}
  				if(r2<0.85){
  					box(lty="solid", col="orange", lwd=2)
  				}
  				if(r2<0.80){
  					box(lty="solid", col="red", lwd=2)
  				}
  				# if this is the last cycle for this fish - save the plot
  				if (i == length(seq_end)){
  					dev.off()
  				}

  			# temp recorded for each measurement cycle
  			temp_mean<-round(mean(d[,temp]),2) # paste this value on the plot
  			temp_max<-round(max(d[,temp]),2)
  			temp_min<-round(min(d[,temp]),2)
  			n_min<-d$time_min[nrow(d)]-d$time_min[1]

  			#oxygen recorded for each measurement cycle
  			O2_mean<-round(mean(d[,Ch], na.rm=TRUE),2)
  			O2_max<-round(max(d[,Ch], na.rm=TRUE),2)
  			O2_min<-round(min(d[,Ch], na.rm=TRUE),2)

  			type<-"SMR"
  			time_frame<-paste("min", round(start,2), "_", round(end,2), sep="")
  			values<-as.data.frame(t(c(time_frame,start, r2, b, m, temp_min, temp_max, temp_mean, O2_min, O2_max, O2_mean, substr(colnames(d[Ch]), start=1, stop=3), DateTime_start,type, n_min)))
  			colnames(values)<-c("time_frame","min_start", "r2" ,"b", "m" , "t_min", "t_max", "t_mean","O2_min", "O2_max", "O2_mean", "Ch", "DateTime_start", "type", "n_min")

  			newdata<-rbind(newdata, values) # add all values of interest to the individual fish specific data file (this is saved in MMR_SMR function)
  			colnames(newdata)<-c("time_frame","min_start", "r2" ,"b", "m" , "t_min", "t_max", "t_mean","O2_min", "O2_max", "O2_mean", "Ch", "DateTime_start", "type", "n_min")
  		}

  	  if(plot_temp==TRUE){
        # run thrugh the loop one more time, but only for temperature
        # no inventory data
        for (i in 1:length(seq_end)){

          start<-seq_st[i]+chop_start # chopping off first x min
    			end<-seq_end[i]-chop_end # chopping off last x min

    			d<-data1[c(which(data1$time_min>start & data1$time_min<end)),] # from teh entire files, get only the data section we need
    			# print(plotname_temperature)
          ### only temp plots:
      			# if this is the first cycle for this fish - start a new plot
      			if (i == 1){

      				if(length(seq_st)<=100){
      					png(plotname_temperature, width=40, height=40, units="in",  res=200)
      					par(mfrow=c(10,10)) # This will fit 100 plots.
      				}
      				if(length(seq_st)>100 & length(seq_st)<=150){
      					png(plotname_temperature, width=40, height=60, units="in",  res=200)
      					par(mfrow=c(15,10)) # This will fit 150 plots
      				}
      				if(length(seq_st)>150 & length(seq_st)<=225){
      					png(plotname_temperature, width=50, height=50, units="in",  res=200)
      					par(mfrow=c(15,15)) # This will fit 150 plots
      				}
          		if(length(seq_st)>=225 & length(seq_st)<=375){
      			  	png(plotname_temperature, width=50, height=70, units="in",  res=200)
      				  par(mfrow=c(15,25)) # This will fit 375 plots
      		    }
      		    if(length(seq_st)>375 ){
      			   	png(plotname_temperature, width=65, height=80, units="in",  res=200)
      			  	par(mfrow=c(18,35)) # This will > 525 plots
      		    }

      			}

      				# plotting each segment with slope, r2, and equations on it
      			plot(d[,temp]~time_min, d=d, ylab="Temperature (C)", xlab="Time (min)", col="darkgreen")
      				# if this is the last cycle for this fish - save the plot
      				if (i == length(seq_end)){
      					dev.off()
      				}
      ### end temp plots
         } # end of the temp plot loop
  	  }

  	}

    # yes inventory data-----
  	if(nrow(inv.data)>0){
  	  if(Ch==4){
  	    cat(paste("Ch1: Data are cleaned for quality sections only \n"))
  	  }else{
  	    cat(paste("Ch",Ch-4 ,": Data are cleaned for quality sections only \n", sep =""))
  	  }

  	  cols = c(3, 4, 5, 6)
  	  for (i in 1:length(cols)){
	       inv.data[, i]<-as.numeric(as.character(inv.data[, i]))
  	  }

  		inv.data.clean<-inv.data[!(inv.data$type=="slope_analysis"),]

  		for (i in 1:length(seq_end)){
  			# if this is the first cycle for this fish - start a new data file
  				# if(as.character(newdata[1,1])=="new"){
  				# 	newdata<-matrix(ncol=15,nrow=0)
  				# 	colnames(newdata)<-c("time_frame","min_start", "r2" ,"b", "m" , "t_min", "t_max", "t_mean","O2_min", "O2_max", "O2_mean", "Ch", "DateTime_start", "type", "n_min")
  				# }

  				# if this is the first cycle for this fish - start a new plot
  				if (i == 1){

  					if(length(seq_st)<=100){
  						png(plotname, width=40, height=40, units="in",  res=200)
  						par(mfrow=c(10,10)) # This will fit 100 plots.
  					}
  					if(length(seq_st)>100 & length(seq_st)<=150){
  						png(plotname, width=40, height=60, units="in",  res=200)
  						par(mfrow=c(15,10)) # This will fit 150 plots
  					}
  					if(length(seq_st)>150){
  						png(plotname, width=50, height=50, units="in",  res=200)
  						par(mfrow=c(15,15)) # This will fit 150 plots
  					}

  				}

  			start<-seq_st[i]+chop_start # chopping off first x min from the beginning of each slopes
  			end<-seq_end[i]-chop_end # chopping off last x min from the end of each slope

  				n2<-which((inv.data.clean[,5]>=start-1 & inv.data.clean[,5]<=end+1) |
  				            (inv.data.clean[,6]>=start-1 & inv.data.clean[,6]<=end+1))

  				if(length(n2)==0){
  					cycle_use<-"use full cycle"

  				}else{

  				  if(inv.data.clean[n2,4] == 1){
  				    cycle_use<-"skip cycle"
  				  }else{
  				    cycle_use<-"use cleaned cycle"
  				  }

  				  if(inv.data.clean[n2,4]==0){
  						inv.data.clean[n2,5]<-start
  					}

  					d_clean<-data1[c(which(data1$time_min>inv.data.clean[n2,5] & data1$time_min<inv.data.clean[n2,6])),]

  				}

  			if (cycle_use=="use full cycle" | cycle_use=="use cleaned cycle" | cycle_use=="skip cycle"){
  				d0<-data1[c(which(data1$time_min>start & data1$time_min<end)),] # from the entire files, get only the data section we need / not cleaned no need to be cleaned

  					if(cycle_use=="use cleaned cycle"){
  						# cleaned section only
  						DateTime_start_clean<-as.character(d_clean$DateTime[1])
  						lm_coef_clean<-coef(lm(d_clean[,Ch]~d_clean$time_min)) # get linear regression fit
  						r2_clean<-round(summary(lm(d_clean[,Ch]~d_clean$time_min))$r.squared,3) # get r2
  						m_clean<-round(lm_coef_clean[2],5) # get slope
  						b_clean<-round(lm_coef_clean[1],2) # get intercept
  						n_min_clean<-d_clean$time_min[nrow(d_clean)]-d_clean$time_min[1]

  						# temp recorded for each measure cycle
  						temp_mean_clean<-round(mean(d_clean[,temp]),2) # paste this value on the plot
  						temp_max_clean<-round(max(d_clean[,temp]),2)
  						temp_min_clean<-round(min(d_clean[,temp]),2)

  						# O2 recorded for each measure cycle
  						O2_mean_clean<-round(mean(d_clean[,Ch]),2)
  						O2_max_clean<-round(max(d_clean[,Ch]),2)
  						O2_min_clean<-round(min(d_clean[,Ch]),2)

  						time_frame<-paste("min", round(inv.data.clean[n2,5],2), "_", round(inv.data.clean[n2,6],2), sep="")
  						type="SMR"

  						if(!cycle_use=="skip cycle"){
  						  values<-as.data.frame(t(c(time_frame, inv.data.clean[n2,5], r2_clean, b_clean, m_clean, temp_min_clean, temp_max_clean, temp_mean_clean,O2_min_clean, O2_max_clean, O2_mean_clean, substr(colnames(d_clean[Ch]), start=1, stop=3), DateTime_start_clean, type, n_min_clean)))
  						  colnames(values)<-c("time_frame","min_start", "r2" ,"b", "m" , "t_min", "t_max", "t_mean","O2_min", "O2_max", "O2_mean", "Ch", "DateTime_start", "type", "n_min")
  						}

  					}else{
  						# full data / not cleaned
  						DateTime_start0<-as.character(d0$DateTime[1])
  						lm_coef0<-coef(lm(d0[,Ch]~d0$time_min)) # get linear regression fit
  						r20<-round(summary(lm(d0[,Ch]~d0$time_min))$r.squared,3) # get r2
  						m0<-round(lm_coef0[2],5) # get slope
  						b0<-round(lm_coef0[1],2) # get intercept
  						n_min0<-d0$time_min[nrow(d0)]-d0$time_min[1]
  						# temp recorded for each measure cycle
  						temp_mean0<-round(mean(d0[,temp]),2) # paste this value on the plot
  						temp_max0<-round(max(d0[,temp]),2)
  						temp_min0<-round(min(d0[,temp]),2)

  						# O2 recorded for each measure cycle
  						O2_mean0<-round(mean(d0[,Ch]),2) # paste this value on the plot
  						O2_max0<-round(max(d0[,Ch]),2)
  						O2_min0<-round(min(d0[,Ch]),2)

  						time_frame<-paste("min", round(start,2), "_", round(end,2), sep="")
  						type="SMR"

    						if(!cycle_use=="skip cycle"){
    						  values<-as.data.frame(t(c(time_frame, start, r20, b0, m0, temp_min0, temp_max0, temp_mean0,O2_min0, O2_max0, O2_mean0, substr(colnames(d0[Ch]), start=1, stop=3), DateTime_start0, type, n_min0)))
    						  colnames(values)<-c("time_frame","min_start", "r2" ,"b", "m" , "t_min", "t_max", "t_mean","O2_min", "O2_max", "O2_mean", "Ch", "DateTime_start", "type", "n_min")
    						}
  # 						  if(cycle_use == "skip cycle" & i == 1){
  #                 # if this is the first cycle for this fish and the newdata is empty
  # 						    # reset newdata to "new
  #           	    newdata<-as.data.frame("new")
  # 						  }
  						}

  				# plotting each segment with slope, r2, and equations on it
  				plot(d0[,Ch]~time_min, d=d0, ylab=expression(paste(O2~(mg~L^{-1}))),xlab="Time (min)")
  					if(cycle_use=="use cleaned cycle"){
  						rect(xleft=inv.data.clean[n2,5], ybottom=min(d0[,Ch]),xright=inv.data.clean[n2,6],ytop=max(d0[,Ch]) ,col= alpha("grey", 0.3), border=NA)
  						points(d_clean[,Ch]~time_min, d=d_clean, ylab=expression(paste(O2~(mg~L^{-1}))),xlab="Time (min)")
  						abline(lm(d_clean[,Ch]~d_clean$time_min), col="red",lwd=2)
  						mtext(bquote(y == .(lm_coef_clean[2])*x + .(lm_coef_clean[1])), adj=1, padj=0, cex=0.8, line=0) # display equation
  						mtext(bquote(italic(R)^2 == .(format(r2_clean, digits = 3))),adj=1, padj=0, cex=0.8, line=1)

  						if(r2_clean<0.9){
  							box(lty="dashed", col="orange", lwd=2)
  						}
  						if(r2_clean<0.85){
  							box(lty="solid", col="orange", lwd=2)
  						}
  						if(r2_clean<0.80){
  							box(lty="solid", col="red", lwd=2)
  						}

  					}else{
  						abline(lm(d0[,Ch]~d0$time_min), col="red",lwd=2)
  						mtext(bquote(y == .(lm_coef0[2])*x + .(lm_coef0[1])), adj=1, padj=0, cex=0.8, line=0) # display equation
  						mtext(bquote(italic(R)^2 == .(format(r20, digits = 3))),adj=1, padj=0, cex=0.8, line=1)

  						if(cycle_use=="skip cycle"){
  						  legend("center","center", c("EXCLUDING \n from dataset"),  pch=" ", cex=2)
  						  box(lty="solid", col="red", lwd=3)
  						}

  						if(r20<0.9){
  							box(lty="dashed", col="orange", lwd=2)
  						}
  						if(r20<0.85){
  							box(lty="solid", col="orange", lwd=2)
  						}
  						if(r20<0.80){
  							box(lty="solid", col="red", lwd=2)
  						}
  					}

  					# if this is the last cycle for this fish - save the plot
  					if (i == length(seq_end)){
  						dev.off()
  					}

  				if(!cycle_use=="skip cycle"){
    				newdata<-rbind(newdata, values) # add all values of interest to the individual fish specific data file (this is saved in MMR_SMR function)
    				colnames(newdata)<-c("time_frame","min_start", "r2" ,"b", "m" ,
    				                     "t_min", "t_max", "t_mean","O2_min", "O2_max",
    				                     "O2_mean", "Ch", "DateTime_start", "type", "n_min")
  				}
  			}

  		}

  	  # run through the loop one more time, but only for temperature
      # WITH inventory data
  		if(plot_temp==TRUE){
        for (i in 1:length(seq_end)){

        start<-seq_st[i]+chop_start # chopping off first x min from the beginning of each slopes
  			end<-seq_end[i]-chop_end # chopping off last x min from the end of each slope

  			n2<-which((inv.data.clean[,4]>=start-1 & inv.data.clean[,4]<end+1) |
  			            (inv.data.clean[,5]>=start-1 & inv.data.clean[,5]<end+1) |
  			            (inv.data.clean[,6]>=start-1 & inv.data.clean[,6]<end+1))
  				if(length(n2)==0){
  					cycle_use<-"use full cycle"
  				  }else{

  					# replace zeros in the inventory data to real values; in inv. data clean - this is where in invenotory file I added 0 when it just starts from the "start" of teh cycle and ends at the "end" of the cycle
  					if(inv.data.clean[n2,4]==0){
  						inv.data.clean[n2,4]<-start
  					}
  					if(inv.data.clean[n2,5]==0){
  						inv.data.clean[n2,5]<-start
  					}
  					if(inv.data.clean[n2,6]==0){
  						inv.data.clean[n2,6]<-end
  					}
  					if(start==inv.data.clean[n2,5] & end==inv.data.clean[n2,6]){
  						cycle_use<-"skip cycle"
  					}else{
  						cycle_use<-"use cleaned cycle"
  					}

  					d_clean<-data1[c(which(data1$time_min>inv.data.clean[n2,5] & data1$time_min<inv.data.clean[n2,6])),]

  			  } # finalizing start and end values

  			if (cycle_use=="use full cycle" | cycle_use=="use cleaned cycle"){	# this is needed to NOT plot cycles that are discarded
  				d0<-data1[c(which(data1$time_min>start & data1$time_min<end)),] # use the entire files to get only the data section sectin secified by cleaning

      		# if this is the first cycle for this fish - start a new temperature plot
      		if (i == 1){
      			if(length(seq_st)<=100){
      				png(plotname_temperature, width=40, height=40, units="in",  res=200)
      				par(mfrow=c(10,10)) # This will fit 100 plots.
      			}
      			if(length(seq_st)>100 & length(seq_st)<=150){
      				png(plotname_temperature, width=40, height=60, units="in",  res=200)
      				par(mfrow=c(15,10)) # This will fit 150 plots
      			}
      			if(length(seq_st)>150){
      				png(plotname_temperature, width=50, height=50, units="in",  res=200)
      				par(mfrow=c(15,15)) # This will fit 150 plots
      			}
      		}

    			# plotting each segment with slope, r2, and equations on it
    		  plot(d0[,temp]~time_min, d=d0, ylab="Temperature (C)", xlab="Time (min)", col="darkgreen")
    		  # data plotted, not make the grey box for the cleaned section
    			if(cycle_use=="use cleaned cycle"){
    				rect(xleft=inv.data.clean[n2,5], ybottom=min(d0[,temp]),xright=inv.data.clean[n2,6],ytop=max(d0[,temp]) ,col= alpha("grey", 0.3), border=NA)
    			}
    			# if this is the last cycle for this fish - save the plot
    		  if (i == length(seq_end)){
    				dev.off()
    			}
  			}

      }
  		}

  	}

  	return(newdata)
    # return(N_Ch)

  }

  	newdata<-matrix(ncol=15,nrow=0)
  	colnames(newdata)<-c("time_frame","min_start", "r2" ,"b", "m" , "t_min", "t_max", "t_mean","O2_min", "O2_max", "O2_mean", "Ch", "DateTime_start", "type", "n_min")

		if(file.exists(data) | file.exists(paste("./csv_files/", data, sep=""))){ # after running through RMRrepeat - this will be saved in csv input files
  	  if(file.exists(paste("./csv_files/", data, sep=""))){
        data1 <- read.csv(paste("./csv_files/", data, sep=""))
        data1 <- data1[complete.cases(data1[ , c('time_sec')]), ]
      }
      if(file.exists(data)){
        data1<-read.csv(data)
      }
  	}else{
      stop_function<-TRUE
      if(stop_function){
        stop("Cannot locate the indicated data data file.")
      }
    }

	# current directory if no specific folders are used
	if(local_path){
	  plotname1<-paste(gsub('.{4}$', '', data), "_all_ChO2.png", sep='')

	  plotname2<-paste( gsub('.{4}$', '', data), "_full.png", sep='')
  	plotname2.1<-paste( gsub('.{4}$', '', data), "_Ch1.png", sep='')
  	plotname2.2<-paste( gsub('.{4}$', '', data), "_Ch2.png", sep='')
  	plotname2.3<-paste( gsub('.{4}$', '', data), "_Ch3.png", sep='')
  	plotname2.4<-paste( gsub('.{4}$', '', data), "_Ch4.png", sep='')

  	plotname_full<-paste(gsub('.{4}$', '', data), "ALL_TOGETHER.png", sep='')

  	filename<-paste(gsub('.{4}$', '', data), "_analyzed.csv", sep='') # save file in the current SMR wd
	  filename2<-paste(gsub('.{4}$', '', data), "_analyzed.csv", sep='') # save in the EPOC_AS etc folder for final SMR and fiull EPOC analysis

	}

	if(!local_path | background_data | sda_data){

	  if(background_data){
      plotname1<-paste("./BACK_RESP/plots_summary_respo/", gsub('.{4}$', '', data), "_all_ChO2.png", sep='')
  	  plotname2<-paste("./BACK_RESP/plots_channel/", gsub('.{4}$', '', data), "_full.png", sep='')

  	  plotname2.1<-paste("./BACK_RESP/plots_channel/", gsub('.{4}$', '', data), "_Ch1.png", sep='')
    	plotname2.2<-paste("./BACK_RESP/plots_channel/", gsub('.{4}$', '', data), "_Ch2.png", sep='')
    	plotname2.3<-paste("./BACK_RESP/plots_channel/", gsub('.{4}$', '', data), "_Ch3.png", sep='')
    	plotname2.4<-paste("./BACK_RESP/plots_channel/", gsub('.{4}$', '', data), "_Ch4.png", sep='')

    	plotname_full<-paste("./BACK_RESP/plots_all_together/",gsub('.{4}$', '', data), "ALL_TOGETHER.png", sep='')
    	filename<-paste("./BACK_RESP/csv_analyzed/", gsub('.{4}$', '', data), "_analyzed.csv", sep='')
        if(!sda_data){
      	  filename2<-paste("./MMR_SMR_AS_EPOC/csv_input_files/", gsub('.{4}$', '', data), "_analyzed.csv", sep='') # save in the EPOC_AS etc folder for final SMR and full EPOC analysis
        }else{
      	  filename2<-paste("./SDA/csv_input_files/", gsub('.{4}$', '', data), "_analyzed.csv", sep='') # save in the EPOC_AS etc folder for final SMR and full EPOC analysis
        }

	  }else{

      plotname1<-paste("./SMR/plots_summary_respo/", gsub('.{4}$', '', data), "_all_ChO2.png", sep='')
  	  plotname2<-paste("./SMR/plots_channel/", gsub('.{4}$', '', data), "_full.png", sep='')

  	  plotname2.1<-paste("./SMR/plots_channel/", gsub('.{4}$', '', data), "_Ch1.png", sep='')
    	plotname2.2<-paste("./SMR/plots_channel/", gsub('.{4}$', '', data), "_Ch2.png", sep='')
    	plotname2.3<-paste("./SMR/plots_channel/", gsub('.{4}$', '', data), "_Ch3.png", sep='')
    	plotname2.4<-paste("./SMR/plots_channel/", gsub('.{4}$', '', data), "_Ch4.png", sep='')

    	plotname_full<-paste("./SMR/plots_all_together/",gsub('.{4}$', '', data), "ALL_TOGETHER.png", sep='')
    	filename<-paste("./SMR/csv_analyzed/", gsub('.{4}$', '', data), "_analyzed.csv", sep='') # save file in the current SMR wd
        if(!sda_data){
      	  filename2<-paste("./MMR_SMR_AS_EPOC/csv_input_files/", gsub('.{4}$', '', data), "_analyzed.csv", sep='') # save in the EPOC_AS etc folder for final SMR and full EPOC analysis      }
          }else{
      	  filename2<-paste("./SDA/csv_input_files/", gsub('.{4}$', '', data), "_analyzed.csv", sep='') # save in the EPOC_AS etc folder for final SMR and full EPOC analysis      }
        }
	  }

	}

	if(plot_temp){
	  if(local_path | !dir.exists("SMR")){
      if(!dir.exists("plots_channel_temperature")){
        dir.create(file.path("./plots_channel_temperature"), recursive = TRUE)
      }
	  }else{
      dir.create(file.path("./SMR/","plots_channel_temperature"), recursive = TRUE)
	  }
    message("Creating new directory: find the return files in subfolder: \"/plots_channel_temperature")
	}


	if(is.null(inventory_data)){
		data2<-as.data.frame(matrix(nrow=0, ncol=0))
	}else{

	   if(file.exists(inventory_data) | file.exists(paste("./SMR/csv_files/", inventory_data, sep=""))){ # after running through RMRrepeat - this will be saved in csv input files
  	  if(file.exists(paste("./SMR/csv_files/", inventory_data, sep=""))){
        data2<-read.csv(paste("./SMR/csv_files/", inventory_data, sep=""))
      }
      if(file.exists(inventory_data)){
        data2<-read.csv(inventory_data)
      }
  	}else{
      stop_function<-TRUE
      if(stop_function){
        stop("Cannot locate the indicated inventory_data data file.")
      }
    }

		 if (plot_temp ==TRUE){
        message("Temperature plots are not saved when using inventory files / cleaning data")
		 }
	}

	if(as.character(data1$Ch1_O2[1])=="--- " || as.character(data1$Ch1_O2[1])=="---" || is.na(data1$Ch1_O2[1])){
		data1$Ch1_O2<-0
	}
	if(as.character(data1$Ch2_O2[1])=="--- " || as.character(data1$Ch2_O2[1])=="---" || is.na(data1$Ch2_O2[1])){
		data1$Ch2_O2<-0
	}
	if(as.character(data1$Ch3_O2[1])=="--- " || as.character(data1$Ch3_O2[1])=="---" || is.na(data1$Ch3_O2[1])){
		data1$Ch3_O2<-0
	}
	if(as.character(data1$Ch4_O2[1])=="--- " || as.character(data1$Ch4_O2[1])=="---" || is.na(data1$Ch4_O2[1])){
		data1$Ch4_O2<-0
	}

	data1$Ch1_temp<-as.numeric(as.character(data1$Ch1_temp))
	data1$Ch1_O2<-as.numeric(as.character(data1$Ch1_O2))
	data1$Ch2_O2<-as.numeric(as.character(data1$Ch2_O2))
	data1$Ch3_O2<-as.numeric(as.character(data1$Ch3_O2))
	data1$Ch4_O2<-as.numeric(as.character(data1$Ch4_O2))

	if (N_Ch==8){
  	data1$Ch2_temp<-as.numeric(as.character(data1$Ch2_temp))
  	data1$Ch3_temp<-as.numeric(as.character(data1$Ch3_temp))
  	data1$Ch4_temp<-as.numeric(as.character(data1$Ch4_temp))
  		#
  }

	#oxygen
	O2_mean1<-mean(data1$Ch1_O2, na.rm=TRUE)
	O2_max1<-max(data1$Ch1_O2, na.rm=TRUE)
	O2_min1<-min(data1$Ch1_O2, na.rm=TRUE)

	O2_mean2<-mean(data1$Ch2_O2, na.rm=TRUE)
	O2_max2<-max(data1$Ch2_O2, na.rm=TRUE)
	O2_min2<-min(data1$Ch2_O2, na.rm=TRUE)

	O2_mean3<-mean(data1$Ch3_O2, na.rm=TRUE)
	O2_max3<-max(data1$Ch3_O2, na.rm=TRUE)
	O2_min3<-min(data1$Ch3_O2, na.rm=TRUE)

	O2_mean4<-mean(data1$Ch4_O2, na.rm=TRUE)
	O2_max4<-max(data1$Ch4_O2, na.rm=TRUE)
	O2_min4<-min(data1$Ch4_O2, na.rm=TRUE)

	message (paste("Minimal O2: ", "[Ch1: ",
	               round(O2_min1,2),"] [Ch2: ",
	               round(O2_min2,2), "] [Ch3: ",
	               round(O2_min3,2), "] [Ch4: ",
	               round(O2_min4,2), "] (mgO2/L)", sep=""))

	data1$date<-as.character(data1$date)
  # herehere 24 h clock
	data1$time<-as.character(data1$time) # this
	print(data1$time)

	DateTime<-strptime(paste(data1$date, data1$time), format = date_format[1], tz = date_format[2])

  if(is.na(DateTime[1])){
    message(paste("DateTime is NA, likely wrong format provided: is it", date_format[1], "?"))
  }

	data1$DateTime<-DateTime
	data1$hr<-round(data1$time_sec/3600,2)
	data1$hr2<-round(data1$time_sec/3600,0)
	data1$time_min<-round(data1$time_sec/60,2)

	png(plotname1, width=25, height=18,units="in",  res=200)
	par(mfrow=c(5,1))
	# cannel 1
	plot(data1$time_min, data1$Ch1_O2, main=paste(data1$DateTime[1], "Channel1"))
	abline(h=O2_min1, lty=2, col="darkred")
	# abline(h=hr_min_mean, lty=1, col="darkred")
	# mtext(paste("temp_min=", temp_min, "temp_max=", temp_max),adj=1, padj=0, cex=0.8, line=1)

	# channel 2
	plot(data1$time_min, data1$Ch2_O2, main=paste(data1$DateTime[1], "Channel2"))
	abline(h=O2_min2, lty=2, col="darkred")
	# abline(h=hr_min_mean, lty=1, col="darkred")

	# channel 3
	plot(data1$time_min, data1$Ch3_O2, main=paste(data1$DateTime[1], "Channel3"))
	abline(h=O2_min3, lty=2, col="darkred")
	# abline(h=hr_min_mean, lty=1, col="darkred")

	# channel 4
	plot(data1$time_min, data1$Ch4_O2, main=paste(data1$DateTime[1], "Channel4"))
	abline(h=O2_min4, lty=2, col="darkred")
	# abline(h=hr_min_mean, lty=1, col="darkred")

	# temperature trace
	plot(data1$time_min, data1$Ch1_temp, main= "Temperature", col="grey", ylim=c(5, 35))
	if (N_Ch>4){
	  	points(data1$time_min, data1$Ch2_temp, main= "Temperatures", col="brown")
	  	points(data1$time_min, data1$Ch3_temp, main= "Temperature", col="orange")
    	points(data1$time_min, data1$Ch4_temp, main= "Temperature", col="purple")
	  }
	dev.off()

	# make this interval adjustable
  t_m<-max(data1$time_min, na.rm = T) # the max min in the file

	if(first_cycle == "flush"){
	    seq_st<-seq(cycle_start, t_m, cycle_end)
	    seq_end<-seq(cycle_end, t_m, cycle_end)
	  } else if(first_cycle == "measure"){
	    seq_st<-seq(0, t_m, cycle_end)
	    seq_end<-seq(cycle_end-cycle_start, t_m, cycle_end)
	  } else {
      stop_function<-TRUE
      if(stop_function) {
        stop("Settings missing: must provide whether the file starts with flush or measurement cycle")
      }
	  }


	rows<-c(4,6,7,8) # the order Ch1, Ch2, Ch3, Ch4

	if(N_Ch==8){
  	rows_temp<-c(5,9,10,11) # the order Ch1, Ch2, Ch3, Ch4
	  } else {
	 	rows_temp<-c(5,5,5,5) # the order Ch1, Ch2, Ch3, Ch4
	}

		if (!data1$Ch1_O2[1]==0){

		  if(!nrow(data2)==0 ){
  			inv.data<-data2[which(grepl(x = data, pattern = as.character(data2[,2]))
  			                      & as.numeric(data2[,3])==1),]
		  }else{
        inv.data <- data2
      }

			newdata<-Channel(Ch=rows[1], temp=rows_temp[1], seq_st, seq_end,
				                 plotname2.1, data1, chop_start, chop_end,
				                 inv.data = inv.data, newdata, plot_temp, N_Ch)

		}

		if (!data1$Ch2_O2[1]==0){

		  if(!nrow(data2)==0 ){
  			inv.data<-data2[which(grepl(x = data, pattern = as.character(data2[,2]))
  			                      & as.numeric(data2[,3])==2),]
		  }else{
        inv.data<-data2
		  }

			newdata<-Channel(Ch=rows[2], temp=rows_temp[2], seq_st, seq_end,
				                 plotname2.1, data1, chop_start, chop_end,
				                 inv.data = inv.data, newdata, plot_temp, N_Ch)

		}

		if (!data1$Ch3_O2[1]==0){

		  if(!nrow(data2)==0 ){
  			inv.data<-data2[which(grepl(x = data, pattern = as.character(data2[,2]))
  			                      & as.numeric(data2[,3])==3),]
		  }else{
        inv.data<-data2
		  }

			newdata<-Channel(Ch=rows[3], temp=rows_temp[3], seq_st, seq_end,
				                 plotname2.1, data1, chop_start, chop_end,
				                 inv.data = inv.data, newdata, plot_temp, N_Ch)

		}
		if (!data1$Ch4_O2[1]==0){

		  if(!nrow(data2)==0 ){
  			inv.data<-data2[which(grepl(x = data, pattern = as.character(data2[,2]))
  			                      & as.numeric(data2[,3])==4),]
		  }else{
        inv.data<-data2
		  }

			newdata<-Channel(Ch=rows[4], temp=rows_temp[4], seq_st, seq_end,
				                 plotname2.1, data1, chop_start, chop_end,
				                 inv.data = inv.data, newdata, plot_temp, N_Ch)
		}

	newdata$min_start<-as.numeric(as.character(newdata$min_start))
	newdata$m<-abs(as.numeric(as.character(newdata$m)))
	newdata$ID_code<-substr(data,start=6, stop=14) # modify this to your needs

	graphics.off()


	plot<-ggplot(data=newdata, aes(min_start, m, col=Ch) )+
		geom_point()+
		geom_line()+
		theme_classic()+
		ylab("Slope (in absolute value; mg/O2/min)")+
		xlab("Time experiment (min)")+
		theme(legend.title=element_blank())


	png(plotname_full, width=9, height=4, res=200, units='in')
		print(plot)
	dev.off()

	write.csv(file=filename, newdata, row.names=FALSE)
	write.csv(file=filename2, newdata, row.names=FALSE)

}

