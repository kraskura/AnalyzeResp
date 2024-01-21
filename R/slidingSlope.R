

#' Title
#'
#' @param data Dataframe with each measurement present (O2 observed each second )
#' @param Ch Channel that needs to be analyzed (use format: c(1, 2), c(3), etc. up to 4 channels)
#' @param local_path Logical. If TRUE (no default) all returned files will be saved in the local working directory. Can also provide a path if this function is run independently
#' @param N_Ch Number of channels for the oxygen meter
#' @return The output from \code{\link{print}}
#' @export
#'
#'
slidingSlope<-function(data,
                       Ch,
                       local_path,
                       N_Ch){

    d<-read.csv(data)
    d$time_min<-round(d$time_sec/60,2)

    #  binding global variables locally to the function.
    s60_1<-s90_1<-s120_1<-s180_1<-cycle_mmr<-m<-r2-NULL

		# for the sliding ones: cycle_type=MMR_slide
		# for the sliding ones: cycle_start=the min on when the slide starts
		# for the sliding ones: cycle_end=the min on when the slide ends
		# for the sliding ones: delay_min=NA
		# for the sliding ones: cycle_mmr= 60,90,120,180 (indicate the sliding length

  # cycle through n channels

  for(n in 1:length(Ch)){
		if(n == 1){
      newdata_mmr<-matrix(ncol = 11, nrow =0)
      colnames(newdata_mmr)<-c("cycle_type","cycle_start","cycle_end",
                              "cycle_mmr", "r2" ,"m", "b" , "t_min", "t_max",
                              "t_mean", "Ch")
		}

    if(Ch[n] == 1){
      # channel 1, column 4
      r <- 4    # temperature only
      r_temp<-5
    }
    if(Ch[n] == 2){
      # channel 1, column 4
      r <- 6
      if(N_Ch==8){
      	r_temp<-9 # the order Ch1, Ch2, Ch3, Ch4
    	  }else {
    	 	r_temp<-5 # the order Ch1, Ch2, Ch3, Ch4
    	}
    }
    if(Ch[n] == 3){
      # channel 1, column 4
      r <- 7
      if(N_Ch==8){
      	r_temp<-10 # the order Ch1, Ch2, Ch3, Ch4
    	  }else {
    	 	r_temp<-5 # the order Ch1, Ch2, Ch3, Ch4
    	}
    }
    if(Ch[n] == 4){
      # channel 1, column 4
      r = 8
      if(N_Ch==8){
      	r_temp<-11 # the order Ch1, Ch2, Ch3, Ch4
    	  }else {
    	 	r_temp<-5 # the order Ch1, Ch2, Ch3, Ch4
    	}
    }

    newdata_set<-matrix(ncol=11, nrow=0)
    colnames(newdata_set)<-c("cycle_type", "cycle_start","cycle_end",
                             "cycle_mmr", "r2" ,"m", "b" , "t_min", "t_max",
                             "t_mean", "Ch")

    newdata_setMean<-matrix(ncol=14, nrow=0)
    colnames(newdata_setMean)<-c("cycle_type",
                                 "cycle_startMean","cycle_endMean", "cycle_mmrMean",
                                 "min_m", "mean_m" ,"sd_m",
                                 "min_r2", "mean_r2" ,"sd_r2",
                                 "t_min", "t_max", "t_mean",
                                 "Ch")

			cycle_type="MMR_slide"

			for (i in 2:nrow(d)){
			  d$time_sec[i]<-d$time_sec[i]-d$time_sec[1]
			}
			d$time_sec[1]<-0

			endtime<-d$time_sec[nrow(d)]
			starttime<-d$time_sec[1]

			if(endtime>180){
				length_slide<-c(60,90,120,180)

				for (k in 1:length(length_slide)){

					string_eval <- sprintf("
						for (i in seq(starttime,(endtime-%s), by=%s)){

							l<-i+%s

							data_set<-d[d$time_sec>i & d$time_sec<l, ]
							DateTime_start<-as.character(data_set$DateTime[1])

							fit_set<-lm(data_set[,r]~data_set$time_min)
							lm_coef_set <- round(coef(fit_set), 5)
							r2_set<-round(as.numeric(as.character(summary(fit_set)$adj.r.squared)),2)

							interc_set<-round(lm_coef_set[1],2)
							slope_set<-round(lm_coef_set[2],6)

							temp_mean_set<-round(mean(data_set[,r_temp]),2)
							temp_max_set<-round(max(data_set[,r_temp]),2)
							temp_min_set<-round(min(data_set[,r_temp]),2)

							cycle_start<-data_set$time_min[1]
							cycle_end<-data_set$time_min[nrow(data_set)]
							cycle_mmr<-%s

							values_set<-as.data.frame(t(c(cycle_type, cycle_start, cycle_end,  cycle_mmr, r2_set ,slope_set, interc_set , temp_min_set, temp_max_set, temp_mean_set, Ch[n])))
							colnames(values_set)<-c(\"cycle_type\", \"cycle_start\",\"cycle_end\", \"cycle_mmr\", \"r2\" ,\"m\", \"b\" , \"t_min\", \"t_max\", \"t_mean\", \"Ch\")

							newdata_set<-rbind(newdata_set, values_set)

						}


					", length_slide[k], 1, length_slide[k],length_slide[k] )
					eval(parse(text = string_eval))
				}


				for (i in 1:length(length_slide)){

					subsetName<-paste("s",length_slide[i],"_1", sep="")
					subsetD<-newdata_set[newdata_set$cycle_mmr ==length_slide[i],]

					assign(subsetName, subsetD)

				}

				list_slideSets<-list(s60_1, s90_1, s120_1, s180_1)


    		newdata_set$cycle_mmr<- as.numeric(as.character(newdata_set$cycle_mmr))

    		if(local_path){
    		  filename_set<-paste( gsub('.{4}$', '', data), "_SLIDINGset", Ch[n] , ".csv", sep='')
    		  filename_setMean<-paste( gsub('.{4}$', '', data), "_SLIDINGsetMean", Ch[n] , ".csv", sep='')
    		  plotnameMean<-paste( gsub('.{4}$', '', data), "_SLIDING_MMRanalysis", Ch[n] , ".png", sep='')
    		  plotname<-paste( gsub('.{4}$', '', data), "_slidingSlope", Ch[n] , ".png", sep='')

    		}else{
    		  if(!is.character(local_path) & dir.exists("MMR")){
            filename_set<-paste("./MMR/Ch[n]annel_sliding_sets/", gsub('.{4}$', '', data), "_SLIDINGset", Ch[n] , ".csv", sep='')
      		  filename_setMean<-paste("./MMR/Channel_sliding_sets/", gsub('.{4}$', '', data), "_SLIDINGsetMean", Ch[n] , ".csv", sep='')
      		  plotnameMean<-paste("./MMR/Channel_plots_MMRanalysis/", gsub('.{4}$', '', data), "_SLIDING_MMRanalysis", Ch[n] , ".png", sep='')
      		  plotname<-paste("./MMR/Channel_plots_MMRanalysis/", gsub('.{4}$', '', data), "_slidingSlope", Ch[n] , ".png", sep='')
    		  }
          if(is.character(local_path)){
            message(paste("MMR sliding window analysis: returns are saved at a specified working directory: \n", local_path))

            filename_set<-paste(local_path, "/Channel_sliding_sets/", gsub('.{4}$', '', data), "_SLIDINGset", Ch[n] , ".csv", sep='')
      		  filename_setMean<-paste(local_path, "/Channel_sliding_sets/", gsub('.{4}$', '', data), "_SLIDINGsetMean", Ch[n] , ".csv", sep='')
      		  plotnameMean<-paste(local_path, "/Channel_plots_MMRanalysis/", gsub('.{4}$', '', data), "_SLIDING_MMRanalysis", Ch[n] , ".png", sep='')
      		  plotname<-paste(local_path, "/Channel_plots_MMRanalysis/", gsub('.{4}$', '', data), "_slidingSlope", Ch[n] , ".png", sep='')
          }
        }

    		write.csv(newdata_set, file=filename_set, row.names=FALSE)

    		newdata_set<-as.data.frame(newdata_set)
    		newdata_set$cycle_mmr<-factor(newdata_set$cycle_mmr)
    		newdata_set$m<-as.numeric(as.character(newdata_set$m))
    		newdata_set$r2<-as.numeric(as.character(newdata_set$r2))

    		yScaleAdjustVal<-1+abs(mean((newdata_set$m)))

    		png(filename = plotnameMean, height=5, width=6, res=200, units="in")
          print(ggplot(data=newdata_set, aes(x=cycle_mmr, y=m))+
          theme_light()+
          # geom_hline(yintercept = -1*(abs(mean((newdata_set$m)))), color = "darkred", alpha=0.6, lty=2)+
          xlab("Time of sliding windows (sec)")+
          scale_y_continuous(name = "Oxygen decrease slope - proxy MMR", sec.axis = sec_axis(trans = ~.+yScaleAdjustVal, breaks = c(1, 0.99, 0.95, 0.8, 0.5), name = expression(R^2)))+
    		  geom_boxplot(width = 0.5, mapping = aes(x=cycle_mmr, y=r2-yScaleAdjustVal), fill="white", color = "darkred", position = position_nudge(x = +0.5))+
    		  geom_point( mapping = aes(x=cycle_mmr, y=r2-yScaleAdjustVal), fill="white", color = "darkred", alpha=0.5, position = position_nudge(x = +0.5), pch=1, size=1)+
          geom_boxplot(width = 0.5, fill="white", color = "grey")+
          geom_point(pch=1, size=1)+
          ggtitle(paste(Ch, " - MMR metrics", sep="")))
        dev.off()

    		for (i in 1:length(list_slideSets)){

    			df_s<-list_slideSets[[i]]

    			df_s$r2<-as.numeric(as.character(df_s$r2))
    			df_s$m<-as.numeric(as.character(df_s$m))
    			r2_set_max<-(which(df_s$r2 == max(df_s$r2)))
    			slope_set_max<-(which(df_s$m == min(df_s$m)))# slope is negative steepest slope is the min(slope)

    			if(length(slope_set_max)==1){
    				mmrMax<-df_s[slope_set_max,]
    		  	}else{
    				# selecting the one with highest r2
    				message(paste(df_s$Ch[1], ": Several MMR slopes: n = ", length(slope_set_max), sep=""))
    				maxset<-df_s[c(slope_set_max),]
    					if(var(maxset$r2)==0){
    					message(df_s$Ch[1], ": Steepest slopes have equivalently good R2, use the first slope recorded", sep="")
    						mmrMax<-maxset[1,]
    						# if the same them I am selecting the first one - that is also plotted

    					}else{
    						message(df_s$Ch[1], ": Steepest slopes with best R2 selected", sep="")
    						bestr2<-which(maxset$r2 == max(maxset$r2))
    						mmrMax<-maxset[bestr2,]
    					}
    			}

    			cols<-c(2:6,8:10)

    			newdata_mmr[,cols] <- lapply(newdata_mmr[,cols], as.character)
    			mmrMax[,cols] <- lapply(mmrMax[,cols], as.character)
    			newdata_mmr[,cols] <- lapply(newdata_mmr[,cols], as.numeric)
    			mmrMax[,cols] <- lapply(mmrMax[,cols], as.numeric)

    			newdata_mmr<-rbind(newdata_mmr,mmrMax)
    		}

        newdata_mmr60<-newdata_mmr[newdata_mmr$cycle_mmr == 60, ]
        newdata_mmr90<-newdata_mmr[newdata_mmr$cycle_mmr == 90, ]
        newdata_mmr120<-newdata_mmr[newdata_mmr$cycle_mmr == 120, ]

    		png(filename = plotname, height=5, width=6, res=200, units="in")
				  plot(d[,r]~time_min, data=d, col="grey",
				       ylab=expression(paste(O2~(mg~L^{-1}))), xlab="Time (min)",
				       main="steepest slope analysis");
				  abline(lm(d[,r]~d$time_min), col="black",lwd=2);

				    plotrix::ablineclip(lm(d[which(d$time_min > newdata_mmr60$cycle_start &
				                            d$time_min < newdata_mmr60$cycle_end),r]~
				                    d[which(d$time_min > newdata_mmr60$cycle_start &
				                              d$time_min < newdata_mmr60$cycle_end),"time_min"]), # purple
				               col=rgb(red = 128, green = 0, blue = 128, alpha = 230, maxColorValue = 255),
				               lty = "solid", lwd=4,lend="round",
				               x1=newdata_mmr60$cycle_start, x2=newdata_mmr60$cycle_end);
				    plotrix::ablineclip(lm(d[which(d$time_min > newdata_mmr90$cycle_start &
				                            d$time_min < newdata_mmr90$cycle_end),r]~
				                    d[which(d$time_min > newdata_mmr90$cycle_start &
				                              d$time_min < newdata_mmr90$cycle_end),"time_min"]), # orange
				               col=rgb(red = 255, green = 127, blue = 80, alpha = 230, maxColorValue = 255),
				               lty = "solid", lwd=4, lend="round",
				               x1=newdata_mmr90$cycle_start, x2=newdata_mmr90$cycle_end);
				    plotrix::ablineclip(lm(d[which(d$time_min > newdata_mmr120$cycle_start &
				                            d$time_min < newdata_mmr120$cycle_end),r]~
				                    d[which(d$time_min > newdata_mmr120$cycle_start &
				                              d$time_min < newdata_mmr120$cycle_end),"time_min"]),
				               col=rgb(red = 0, green = 128, blue = 128, alpha = 230, maxColorValue = 255), # teal
				               lty = "solid", lwd=4, lend="round",
				               x1=newdata_mmr120$cycle_start, x2=newdata_mmr120$cycle_end);
          dev.off()
			}else{
		  	message(paste(Ch, ": file < 3 min: no steepest slope analysis"))
			}
    }

		return(newdata_mmr)

	}
