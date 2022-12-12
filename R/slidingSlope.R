

#' Title
#'
#' @param d Dataframe with each measurement present (O2 observed each second )
#' @param Ch Channel
#' @param data.MMR The name of the export files
#' @param r The column to use for oxygen trace
#' @param r_temp The column to use for temperature trace
#' @param newdata_mmr The ongoing dataframe of mmr that will be row-bind building as export
#' @param local_path Logical. If TRUE (no default) all returned files will be saved in the local working directory. Can also provide a path if this function is run independently
#'
#' @return The output from \code{\link{print}}
#' @export
slidingSlope<-function(d,
                       Ch,
                       data.MMR,
                       r,
                       r_temp,
                       newdata_mmr,
                       local_path){

    #  binding global variables locally to the function.
    s60_1<-s90_1<-s120_1<-s180_1<-cycle_mmr<-m<-r2<-NULL

		# for the sliding ones: cycle_type=MMR_slide
		# for the sliding ones: cycle_start=the min on when the slide starts
		# for the sliding ones: cycle_end=the min on when the slide ends
		# for the sliding ones: delay_min=NA
		# for the sliding ones: cycle_mmr= 60,90,120,180 (indicate the sliding length
    newdata_set<-matrix(ncol=12, nrow=0)
    colnames(newdata_set)<-c("cycle_type", "cycle_start","cycle_end",  "cycle_mmr", "r2" ,"m", "b" , "t_min", "t_max", "t_mean", "Ch", "DateTime_start")

    newdata_setMean<-matrix(ncol=14, nrow=0)
    colnames(newdata_setMean)<-c("cycle_type", "cycle_startMean","cycle_endMean", "cycle_mmrMean","min_m", "mean_m" ,"sd_m", "min_r2", "mean_r2" ,"sd_r2", "t_min", "t_max", "t_mean", "Ch")
    #~ 		normalize time to have all files start at 0
			cycle_type="MMR_slide"

			for (i in 2:nrow(d)){
			  d$time_sec[i]<-d$time_sec[i]-d$time_sec[1]
			}
			d$time_sec[1]<-0

			endtime<-d$time_sec[nrow(d)]
			starttime<-d$time_sec[1]

			if(endtime>180){
				# message("FILE > 3 min")
##			length_slide<-c(10,20,30,60,90,120,180)
				length_slide<-c(60,90,120,180)
				# timeDiff_slide<-c(1)

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

							values_set<-as.data.frame(t(c(cycle_type, cycle_start, cycle_end,  cycle_mmr, r2_set ,slope_set, interc_set , temp_min_set, temp_max_set, temp_mean_set, Ch, DateTime_start)))
							colnames(values_set)<-c(\"cycle_type\", \"cycle_start\",\"cycle_end\", \"cycle_mmr\", \"r2\" ,\"m\", \"b\" , \"t_min\", \"t_max\", \"t_mean\", \"Ch\", \"DateTime_start\")

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

					list<-list(s60_1, s90_1, s120_1, s180_1)


    		newdata_set$cycle_mmr<- as.numeric(as.character(newdata_set$cycle_mmr))

    		if(local_path){
    		  filename_set<-paste( gsub('.{4}$', '', data.MMR), "_SLIDINGset", Ch , ".csv", sep='')
    		  filename_setMean<-paste( gsub('.{4}$', '', data.MMR), "_SLIDINGsetMean", Ch , ".csv", sep='')
    		  plotnameMean<-paste( gsub('.{4}$', '', data.MMR), "_SLIDING_MMRanalysis", Ch , ".png", sep='')

    		}else{
    		  if(!is.character(local_path) & dir.exists("MANUAL")){
            filename_set<-paste("./MANUAL/channel_sliding_sets/", gsub('.{4}$', '', data.MMR), "_SLIDINGset", Ch , ".csv", sep='')
      		  filename_setMean<-paste("./MANUAL/channel_sliding_sets/", gsub('.{4}$', '', data.MMR), "_SLIDINGsetMean", Ch , ".csv", sep='')
      		  plotnameMean<-paste("./MANUAL/channel_plots_MMRanalysis/", gsub('.{4}$', '', data.MMR), "_SLIDING_MMRanalysis", Ch , ".png", sep='')
    		  }
          if(is.character(local_path)){
            message(paste("MMR sliding window analysis: returns are saved at a specified working directory: \n", local_path))

            filename_set<-paste(local_path, "/channel_sliding_sets/", gsub('.{4}$', '', data.MMR), "_SLIDINGset", Ch , ".csv", sep='')
      		  filename_setMean<-paste(local_path, "/channel_sliding_sets/", gsub('.{4}$', '', data.MMR), "_SLIDINGsetMean", Ch , ".csv", sep='')
      		  plotnameMean<-paste(local_path, "/channel_plots_MMRanalysis/", gsub('.{4}$', '', data.MMR), "_SLIDING_MMRanalysis", Ch , ".png", sep='')
          }
        }

    		write.csv(newdata_set,file=filename_set, row.names=FALSE)

    		newdata_set<-as.data.frame(newdata_set)
    		## Jan 2022:
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

    		for (i in 1:length(list)){

    			df_s<-list[[i]]

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


    			# auto correlation

    			### create a sliding for the mean slopes
    			# if(nrow(df_s) > 60){
    				# for (t in 1:(nrow(df_s)-60)){
    				#
    				# 	k<-t+60
    				#
    				# 	data_setMean<-df_s[t:k, ]
    				#
    				# 	temp_mean_setMean<-mean(as.numeric(as.character(data_setMean$t_mean)))
    				# 	temp_min_setMean<-min(as.numeric(as.character(data_setMean$t_min)))
    				# 	temp_max_setMean<-max(as.numeric(as.character(data_setMean$t_max)))
    				#
    				# 	cycle_startMean<-as.character(data_setMean$cycle_start[1])
    				# 	cycle_endMean<-as.character(data_setMean$cycle_start[nrow(data_setMean)])
    				# 	cycle_mmrMean<-as.character(data_setMean$cycle_mmr[1])
    				#
    				# 	mean_m<-mean(as.numeric(as.character(data_setMean$m)))
    				# 	min_m<-min(as.numeric(as.character(data_setMean$m)))
    				# 	sd_m<-sd(as.numeric(as.character(data_setMean$m)))
    				#
    				# 	mean_r2<-mean(as.numeric(as.character(data_setMean$r2)))
    				# 	min_r2<-min(as.numeric(as.character(data_setMean$r2)))
    				# 	sd_r2<-sd(as.numeric(as.character(data_setMean$r2)))
    				#
    				# 	# get mean and sd for r2
    				#
    				# 	values_setMean<-as.data.frame(t(c("Mean_1minSlopes", cycle_startMean, cycle_endMean, cycle_mmrMean, min_m, mean_m , sd_m, min_r2, mean_r2, sd_r2, temp_min_setMean, temp_max_setMean, temp_mean_setMean, Ch)))
    				# 	colnames(values_setMean)<-c("cycle_type", "cycle_startMean","cycle_endMean", "cycle_mmrMean","min_m", "mean_m" ,"sd_m", "min_r2", "mean_r2" ,"sd_r2", "t_min", "t_max", "t_mean", "Ch")
    				#
    				# 	newdata_setMean<-rbind(newdata_setMean, values_setMean)
    				#
    				# }
    			# }

    	  	}
      # find the lowest value and add it to a dataframe
    		#
    		# listMean<-split(newdata_setMean, newdata_setMean$cycle_mmrMean)
    		#
    		# for (v in 1:length(listMean)){
    		#   df_sMean<-listMean[[v]]
    		#
    		# 	df_sMean$mean_r2<-as.numeric(as.character(df_sMean$mean_r2))
    		# 	df_sMean$min_r2<-as.numeric(as.character(df_sMean$min_r2))
    		# 	df_sMean$mean_m<-as.numeric(as.character(df_sMean$mean_m))
    		# 	r2_set_maxMean<-(which(df_sMean$mean_r2 == max(df_sMean$mean_r2)))
    		# 	slope_set_maxMean<-(which(df_sMean$mean_m == min(df_sMean$mean_m)))# slope is negative steepest slope is the min(slope)
    		#
    		# 	if(length(slope_set_maxMean)==1){
    		# 		mmrMaxMean<-df_sMean[slope_set_maxMean,]
    		#     }else{
    		# 		# selecting the one with the lowest mo2 included, and then with the highest min r2
    		# 		maxsetMean<-df_sMean[c(slope_set_maxMean),]
    		# 		mmrMaxMean<-maxsetMean[which(maxsetMean$min_r2==max(maxsetMean$min_r2))[1],] # if the same them I am selecting the first one, with the highest r2
    		#   }
    		#
    		# 	cols<-c(9,6,7,11)
    		# 	mmrMaxMean[,cols] %<>% lapply(function(x) as.numeric(as.character(x)))
    		#
    		# 	mmrMaxMean_vals<-as.data.frame(t(c("Mean_1minSlopes", mmrMaxMean[,2],
    		# 	                                   as.character(mmrMaxMean[,3]), as.character( mmrMaxMean[,4]),
    		# 	                                   as.character(round(mmrMaxMean[,9],4)), as.character(round(mmrMaxMean[,6], 4)),
    		# 	                                   paste("SD_slope:",as.character(round(mmrMaxMean[,7], 4))), as.character(mmrMaxMean[,11]),
    		# 	                                   as.character(mmrMaxMean[,12]), as.character(mmrMaxMean[,13]), as.character(mmrMaxMean[,14]),
    		# 	                                   as.character(df_s$DateTime_start[1]))))
    		# 	colnames(mmrMaxMean_vals)<-c("cycle_type", "cycle_start","cycle_end",  "cycle_mmr", "r2" ,"m", "b" , "t_min", "t_max", "t_mean", "Ch", "DateTime_start")
    		#
    		# 	newdata_mmr<-rbind(newdata_mmr, mmrMaxMean_vals)
    		# }

    #     cols<-c(2:13)
    # 	  newdata_setMean[,cols] %<>% lapply(function(x) as.numeric(as.character(x)))
    # 		newdata_setMean[,cols] %<>% lapply(function(x) as.numeric(as.character(x)))
    #
    # 		# best slopes using the singly fastest slope analysis
    # 		bestdata<-newdata_mmr[c(newdata_mmr$cycle_type=="MMR_slide" & newdata_mmr$Ch==Ch), ]
    # 		names(bestdata)[names(bestdata) == 'cycle_mmr'] <- 'cycle_mmrMean'
    #
    # 		bestdata$cycle_mmrMean<-as.numeric(as.character(bestdata$cycle_mmrMean))
    # 		# newdata_setMean$cycle_mmrMean<-factor(newdata_setMean$cycle_mmrMean, levels=c('60','90','120','180'))
    #
    # 		newdata_setMean$cycle_startMean+0.5
    #
    #   	 error_plot_m<-ggplot(data=newdata_setMean, aes(y=mean_m, x=cycle_startMean+0.5))+
    #   	    geom_point(pch=21, size=2, fill="black")+
    #   	    geom_point(aes(y=min_m, x=cycle_startMean+0.5), pch=1, size=1)+
    #   	    geom_errorbar(aes(ymin=mean_m+sd_m, ymax=mean_m-sd_m), width=0.1, size=0.2)+
    #   	    theme_light()+
    #   	    geom_abline(data=bestdata, aes(slope=0, intercept=as.numeric(as.character(m))), lty=2, colour="red", size=1)+
    #   	    ggtitle(paste(Ch, " slope", sep=""))+
    #   	    facet_grid(cycle_mmrMean~.)
    #
    #   	 error_plot_r2<-ggplot(data=newdata_setMean, aes(y=mean_r2, x=cycle_startMean+0.5))+
    #   	    geom_point(pch=21, size=2, color="red", fill="red")+
    #   	    geom_point(aes(y=min_r2, x=cycle_startMean+0.5), pch=1, color="red", size=1)+
    #   	    geom_errorbar(aes(ymin=mean_r2+sd_r2, ymax=mean_r2-sd_r2), width=0.1, size=0.2, colour="red")+
    #   	    theme_light()+
    #   	    ggtitle(paste(Ch, " r2", sep=""))+
    #   	    facet_grid(cycle_mmrMean~.)
    #
    #      png(plotnameMean, height=10, width=5, res=300, units="in")
    #   	   grid.arrange(error_plot_m, error_plot_r2, ncol=1, nrow=2)
    #      dev.off()
    #
    #      write.csv(newdata_setMean,file=filename_setMean, row.names=FALSE)
    #

    		#      png(plotnameMean, height=10, width=5, res=300, units="in")
    		#   	   grid.arrange(mmrMeans_plot_m, error_plot_r2, ncol=1, nrow=2)
    		#      dev.off()


			}else{

		  	message(paste(Ch, ": file < 3 min: no steepest slope analysis"))
			# message(c(starttime, endtime))
			}


		return(newdata_mmr)



	}
