#' Title
#'
#' @param data.MMR The name of the raw input file; .csv, a character.
#' @param cycles The number of measurement cycles recording metabolic rate in an animal (a single number value 1 through 4).
#' @param cycle_start A numeric vector of measurement cycle start times (min). The length of this vector has to be equal to the number of cycles provided, e.g. if cycles = 2, then it may be cycle_start = c(0,10))
#' @param cycle_end A numeric vector of measurement cycle end times (min). The same format as cycle_start
#' @param mmr_Ch1 A number indicated in which cycle the animal in oxygen meter Channel 1 is expected to reach its maximum metabolic rate. Must specify to analyze metabolic rate, or enter 0, which will drop the channel from the analysis.
#' @param mmr_Ch2 See mmr_Ch1 argument for description. This applies to Channel 2
#' @param mmr_Ch3 See mmr_Ch1 argument for description. This applies to Channel 3
#' @param mmr_Ch4 See mmr_Ch1 argument for description. This applies to Channel 4
#' @param clean_Ch1 Specific to MMR cycle for each animal (or cycle). Use to set a time window (min) that specifies a quality section on metabolic rate measurement for sliding window analysis. Default will use full length measurement cycle as set by cycle_start and cycle_end arguments.
#' @param clean_Ch2 See mmr_Ch1 argument for description. This applies to Channel 2
#' @param clean_Ch3 See mmr_Ch1 argument for description. This applies to Channel 3
#' @param clean_Ch4 See mmr_Ch1 argument for description. This applies to Channel 4
#' @param N_Ch The number of FireSting channels. Options include 2, 4, and 8. This argument can be ignored if 2-Channel FireSting was used.
#' @param path Specify the working directory for all output files. This argument has two options, i) files can be automatically organized in newly created local directories. Any character may be used.
#' @param date_format The date format used in the original FireSting files. Must specify one of the following: c("m/d/y","d/m/y","y-m-d")
#' @param inv.data The name of an inventory data dataframe that provides dertails in "cleaning"
#'
#' @return The output from \code{\link{print}}
#' @export
#'
#' @importFrom stats lm
#' @importFrom stats coef
#' @importFrom stats var
#' @import graphics
#' @import grDevices
#' @import scales
#' @import chron
#' @import ggplot2
#' @import utils
#' @import stringr
#' @importFrom dplyr group_by
#' @importFrom dplyr summarize
#' @importFrom pryr %<a-%
#' @importFrom plotrix ablineclip
#' @examples
#' \dontrun{
#' MMR(data.MMR = "data.csv", cycles = 1, cycle_start = c(0), cycle_end = c(5))
#' }
MMR<-function(data.MMR,
              cycles,
              cycle_start,
              cycle_end,
              mmr_Ch1=0,
              mmr_Ch2=0,
              mmr_Ch3=0,
              mmr_Ch4=0,
              clean_Ch1=c(0,0), clean_Ch2=c(0,0), clean_Ch3=c(0,0), clean_Ch4=c(0,0),
              N_Ch = 4, path = ".", date_format = "m/d/y",
              inv.data = NA){

  #  binding global variables locally to the function.
  p<-p1<-p2<-p3<-p4<-p_null<-NULL


  Channel_mmr<-function(data.MMR, dataMMR, cycle_start, cycle_end, clean_start_mmr, clean_end_mmr, Ch_list, cycles, j, rows, rows_temp, newdata_mmr, path, inv.data.clean){

		if(!Ch_list[j]==0){
			mmr_start<-clean_start_mmr[j]
			mmr_end<-clean_end_mmr[j]
			r<-rows[j]
			r_temp<-rows_temp[j]
			cycle_mmr<-Ch_list[j]

			d<-dataMMR[c(which(dataMMR$time_min>(mmr_start) & dataMMR$time_min<mmr_end)),]

			lm_coef<-coef(lm(d[,r]~d$time_min))
			r2<-round(summary(lm(d[,r]~d$time_min))$r.squared,3)

			m<-round(lm_coef[2],5)
			b<-round(lm_coef[1],2)

			clean_cycle2<-NA
			clean_cycle3<-NA
			clean_cycle4<-NA

			### is there cleaning?
  		if(nrow(inv.data.clean)>0){
  		   message(paste("Ch",j, ": INVENTORY data, adjusting non-MMR cycle times", sep=""))
  	      # find cycle 2 cleaning
  		   # print(clean_cycle2)
  	      if(any(inv.data.clean$cycle==2)){
  	         clean_cycle2<-"Y"
  	         # print(clean_cycle2)
  	         # if(!inv.data.clean[which(inv.data.clean$cycle==2),"start"] == 0){
  	           start_c2<-inv.data.clean[which(inv.data.clean$cycle==2),"start"] # min start
  	           end_c2<-inv.data.clean[which(inv.data.clean$cycle==2),"end"]
  	         # }
  	         # else{
  	         #   start_c2<-cycle_start[2]
  	         #   end_c2<-cycle_end[2]
  	         # }
  	         # print(start_c2)
  	      }
  		   # find cycle 3 cleaning
  	      if(any(inv.data.clean$cycle==3)){
  	        # print(inv.data.clean)
  	         clean_cycle3<-"Y"
  	         # if(!inv.data.clean[which(inv.data.clean$cycle==3),"start"] == 0){
    	         start_c3<-inv.data.clean[which(inv.data.clean$cycle==3),"start"] # min start
    	         end_c3<-inv.data.clean[which(inv.data.clean$cycle==3),"end"]
#   	         }
#   	         else{
# 	            start_c3<-0
# 	            end_c3<-0
#   	         }

  	      }
  		  # find cycle 4 cleaning
  	      if(any(inv.data.clean$cycle==4)){
  	         # print(inv.data.clean)
  	         clean_cycle4<-"Y"
  	         # if(!inv.data.clean[which(inv.data.clean$cycle==4),"start"] == 0){
    	         start_c4<-inv.data.clean[which(inv.data.clean$cycle==4),"start"] # min start
    	         end_c4<-inv.data.clean[which(inv.data.clean$cycle==4),"end"]
  	         # }
  	         # else{
  	         #   start_c4<-cycle_start[4]
  	         #   end_c4<-cycle_end[4]
  	         # }
  	      }
  		}

			###  end to cleaning times for cycles

				# print(c(clean_cycle2, start_c2, end_c2, j))

			## these are standard cycles without any cleaning
			if (cycles>=2){
			  if(is.na(clean_cycle2)){
			    # start_c2<-cycle_start[2]
			    # end_c2<-cycle_end[2]
				  d2<-dataMMR[c(which(dataMMR$time_min>cycle_start[2] & dataMMR$time_min<cycle_end[2])),]
			  }else{
			    if(start_c2==0 & end_c2 ==0){
			      d2<-dataMMR[c(which(dataMMR$time_min>cycle_start[2] & dataMMR$time_min<cycle_end[2])),]
			    }else{
			      d2<-dataMMR[c(which(dataMMR$time_min>start_c2 & dataMMR$time_min<end_c2)),]
			    }
			  }
  				lm_coef.2<-coef(lm(d2[,r]~d2$time_min))
  				r2.2<-round(summary(lm(d2[,r]~d2$time_min))$r.squared,3)
  				m.2<-round(lm_coef.2[2],5)
  				b.2<-round(lm_coef.2[1],2)
			}
			if (cycles>=3){
			 if(is.na(clean_cycle3)){
  		   # start_c3<-cycle_start[3]
  		   # end_c3<-cycle_end[3]
			   d3<-dataMMR[c(which(dataMMR$time_min>cycle_start[3] & dataMMR$time_min<cycle_end[3])),]
			 }else{
			   if(start_c3==0 & end_c3 ==0){
			     d3<-dataMMR[c(which(dataMMR$time_min>cycle_start[3] & dataMMR$time_min<cycle_end[3])),]
			   }else{
			     d3<-dataMMR[c(which(dataMMR$time_min>start_c3 & dataMMR$time_min<end_c3)),]
			   }
			 }

				lm_coef.3<-coef(lm(d3[,r]~d3$time_min))
				r2.3<-round(summary(lm(d3[,r]~d3$time_min))$r.squared,3)
				m.3<-round(lm_coef.3[2],5)
				b.3<-round(lm_coef.3[1],2)
			}
			if (cycles>=4){
			  if(is.na(clean_cycle4)){
			    # start_c4<-cycle_start[4]
			    # end_c4<-cycle_end[4]
				  d4<-dataMMR[c(which(dataMMR$time_min>cycle_start[4] & dataMMR$time_min<cycle_end[4])),]
			  }else{
			    if(start_c4==0 & end_c4 ==0){
			    }
			      d4<-as.data.frame(dataMMR[c(which(dataMMR$time_min>start_c4 & dataMMR$time_min<end_c4)),])
			    # print(c(end_c4, start_c4))
			  }
				lm_coef.4<-coef(lm(d4[,r]~d4$time_min))
				r2.4<-round(summary(lm(d4[,r]~d4$time_min))$r.squared,3)
				m.4<-round(lm_coef.4[2],5)
				b.4<-round(lm_coef.4[1],2)
			}

  		if (path == "."){
         plotname<-paste( gsub('.{4}$', '', data.MMR),"_", substr(colnames(d[r]), start=1, stop=3), ".png", sep="")#

      	}else{
      	 plotname<-paste("../channel_plots/", gsub('.{4}$', '', data.MMR),"_", substr(colnames(d[r]), start=1, stop=3), ".png", sep="")#
    	}


			# all MMR file trace
			p %<a-%{
				plot(dataMMR[,r]~time_min, data=dataMMR, ylab=expression(paste(O2~(mg~L^{-1}))),xlab="Time (min)", main="The entire MMR file");
				rect(xleft=mmr_start[1],ybottom=min(dataMMR[,r]),xright=mmr_end[1],ytop=max(dataMMR[,r]) ,col= '#FF003322', border=NA)
			}

			# MMR plots
			if (Ch_list[j]==1 & cycles>=1){

				DateTime_start<-as.character(d$DateTime[1])
				cycle_type<-"MMR"
				Ch<-paste("Ch",j, sep="")

				values_mmr<-as.data.frame(t(c(cycle_type, mmr_start, mmr_end,  cycle_mmr, r2, m, b,
				                              round(min(d[,r_temp]), 4),
				                              round(max(d[,r_temp]), 4),
				                              round(mean(d[,r_temp]), 4),
				                              Ch ,DateTime_start )))
				colnames(values_mmr)<-c("cycle_type","cycle_start","cycle_end",  "cycle_mmr", "r2" ,"m", "b" , "t_min", "t_max", "t_mean", "Ch", "DateTime_start")
				newdata_mmr<-rbind(newdata_mmr, values_mmr)


				print(newdata_mmr)
				newdata_mmr<-slidingSlope(d, Ch, data.MMR, r, r_temp, newdata_mmr, path)
				print(newdata_mmr)


				newdata_mmr60<-newdata_mmr[which(newdata_mmr$Ch == Ch & newdata_mmr$cycle_mmr == 60),]
				newdata_mmr90<-newdata_mmr[which(newdata_mmr$Ch == Ch & newdata_mmr$cycle_mmr == 90),]
				newdata_mmr120<-newdata_mmr[which(newdata_mmr$Ch == Ch & newdata_mmr$cycle_mmr == 120),]

				# print(c(Ch_list[j], Ch))

				p1 %<a-% {
				  plot(d[,r]~time_min, data=d, col="grey", ylab=expression(paste(O2~(mg~L^{-1}))), xlab="Time (min)",main="MMR");
				  abline(lm(d[,r]~d$time_min), col="black",lwd=2);
				  if(nrow(newdata_mmr60)>0){

				    ablineclip(lm(d[which(d$time_min > newdata_mmr60$cycle_start &
				                            d$time_min < newdata_mmr60$cycle_end),r]~
				                    d[which(d$time_min > newdata_mmr60$cycle_start &
				                              d$time_min < newdata_mmr60$cycle_end),"time_min"]), # purple
				               col=rgb(red = 128, green = 0, blue = 128, alpha = 230, maxColorValue = 255),
				               lty = "solid", lwd=4,lend="round",
				               x1=newdata_mmr60$cycle_start, x2=newdata_mmr60$cycle_end);
				    ablineclip(lm(d[which(d$time_min > newdata_mmr90$cycle_start &
				                            d$time_min < newdata_mmr90$cycle_end),r]~
				                    d[which(d$time_min > newdata_mmr90$cycle_start &
				                              d$time_min < newdata_mmr90$cycle_end),"time_min"]), # orange
				               col=rgb(red = 255, green = 127, blue = 80, alpha = 230, maxColorValue = 255),
				               lty = "solid", lwd=4, lend="round",
				               x1=newdata_mmr90$cycle_start, x2=newdata_mmr90$cycle_end);
				    ablineclip(lm(d[which(d$time_min > newdata_mmr120$cycle_start &
				                            d$time_min < newdata_mmr120$cycle_end),r]~
				                    d[which(d$time_min > newdata_mmr120$cycle_start &
				                              d$time_min < newdata_mmr120$cycle_end),"time_min"]),
				               col=rgb(red = 0, green = 128, blue = 128, alpha = 230, maxColorValue = 255), # teal
				               lty = "solid", lwd=4, lend="round",
				               x1=newdata_mmr120$cycle_start, x2=newdata_mmr120$cycle_end);
				  }


				  # points(d[,r]~d$time_min, col="grey");
				  mtext(bquote(y == .(lm_coef[2])*x + .(lm_coef[1])), adj=1, padj=0, cex=0.8, line=0); # display equation
				  mtext(bquote(italic(R)^2 == .(format(r2, digits = 3))),adj=1, padj=0, cex=0.8, line=1)
				}
			}

			if (Ch_list[j]==2 & cycles>=2){

				DateTime_start<-as.character(d$DateTime[1])
				cycle_type<-"MMR"
				Ch<-paste("Ch",j, sep="")
				values_mmr<-as.data.frame(t(c(cycle_type, mmr_start, mmr_end,  cycle_mmr, r2, m, b,
				                              round(min(d[,r_temp]), 4),
				                              round(max(d[,r_temp]), 4),
				                              round(mean(d[,r_temp]), 4),
				                              Ch ,DateTime_start )))
				colnames(values_mmr)<-c("cycle_type","cycle_start","cycle_end", "cycle_mmr", "r2" ,"m", "b" , "t_min", "t_max", "t_mean", "Ch", "DateTime_start")
				newdata_mmr<-rbind(newdata_mmr, values_mmr)

				newdata_mmr<-slidingSlope(d, Ch, data.MMR, r, r_temp, newdata_mmr, path)

				newdata_mmr60<-newdata_mmr[which(newdata_mmr$Ch == Ch & newdata_mmr$cycle_mmr == 60),]
				newdata_mmr90<-newdata_mmr[which(newdata_mmr$Ch == Ch & newdata_mmr$cycle_mmr == 90),]
				newdata_mmr120<-newdata_mmr[which(newdata_mmr$Ch == Ch & newdata_mmr$cycle_mmr == 120),]


				p2 %<a-% {
				  plot(d[,r]~time_min, data=d,  col="grey", ylab=expression(paste(O2~(mg~L^{-1}))),xlab="Time (min)",main="MMR");
				  abline(lm(d[,r]~d$time_min), col="black",lwd=2);

				  if(nrow(newdata_mmr60)>0){
				    ablineclip(lm(d[which(d$time_min > newdata_mmr60$cycle_start &
  				                          d$time_min < newdata_mmr60$cycle_end),r]~
  				                  d[which(d$time_min > newdata_mmr60$cycle_start &
  				                            d$time_min < newdata_mmr60$cycle_end),"time_min"]), # purple
  				             col=rgb(red = 128, green = 0, blue = 128, alpha = 230, maxColorValue = 255),
  				             lty = "solid", lwd=4,lend="round",
  				             x1=newdata_mmr60$cycle_start, x2=newdata_mmr60$cycle_end);
  				  ablineclip(lm(d[which(d$time_min > newdata_mmr90$cycle_start &
  				                          d$time_min < newdata_mmr90$cycle_end),r]~
  				                  d[which(d$time_min > newdata_mmr90$cycle_start &
  				                            d$time_min < newdata_mmr90$cycle_end),"time_min"]), # orange
  				             col=rgb(red = 255, green = 127, blue = 80, alpha = 230, maxColorValue = 255),
  				             lty = "solid", lwd=4, lend="round",
  				             x1=newdata_mmr90$cycle_start, x2=newdata_mmr90$cycle_end);
  				  ablineclip(lm(d[which(d$time_min > newdata_mmr120$cycle_start &
  				                          d$time_min < newdata_mmr120$cycle_end),r]~
  				                  d[which(d$time_min > newdata_mmr120$cycle_start &
  				                            d$time_min < newdata_mmr120$cycle_end),"time_min"]),
  				             col=rgb(red = 0, green = 128, blue = 128, alpha = 230, maxColorValue = 255), # teal
  				             lty = "solid", lwd=4, lend="round",
  				             x1=newdata_mmr120$cycle_start, x2=newdata_mmr120$cycle_end);
				  }

				  # points(d[,r]~d$time_min, col="grey");
				  mtext(bquote(y == .(lm_coef[2])*x + .(lm_coef[1])), adj=1, padj=0, cex=0.8, line=0); # display equation
				  mtext(bquote(italic(R)^2 == .(format(r2, digits = 3))),adj=1, padj=0, cex=0.8, line=1);
				}
			}

			if (Ch_list[j]==3 & cycles>=3){

				DateTime_start<-as.character(d$DateTime[1])
				cycle_type<-"MMR"
				Ch<-paste("Ch",j, sep="")
				values_mmr<-as.data.frame(t(c(cycle_type, mmr_start, mmr_end,  cycle_mmr, r2, m, b,
				                              round(min(d[,r_temp]), 4),
				                              round(max(d[,r_temp]), 4),
				                              round(mean(d[,r_temp]), 4),
				                              Ch ,DateTime_start )))
				colnames(values_mmr)<-c("cycle_type","cycle_start","cycle_end",  "cycle_mmr", "r2" ,"m", "b" , "t_min", "t_max", "t_mean", "Ch", "DateTime_start")
				newdata_mmr<-rbind(newdata_mmr, values_mmr)

				newdata_mmr<-slidingSlope(d, Ch, data.MMR, r, r_temp, newdata_mmr, path)
				newdata_mmr60<-newdata_mmr[which(newdata_mmr$Ch == Ch & newdata_mmr$cycle_mmr == 60),]
				newdata_mmr90<-newdata_mmr[which(newdata_mmr$Ch == Ch & newdata_mmr$cycle_mmr == 90),]
				newdata_mmr120<-newdata_mmr[which(newdata_mmr$Ch == Ch & newdata_mmr$cycle_mmr == 120),]


				p3 %<a-% {
				  plot(d[,r]~time_min, data=d, col="grey", ylab=expression(paste(O2~(mg~L^{-1}))),xlab="Time (min)",main="MMR");
				  abline(lm(d[,r]~d$time_min), col="black",lwd=2);


				  if(nrow(newdata_mmr60)>0){
				    ablineclip(lm(d[which(d$time_min > newdata_mmr60$cycle_start &
  				                          d$time_min < newdata_mmr60$cycle_end),r]~
  				                  d[which(d$time_min > newdata_mmr60$cycle_start &
  				                            d$time_min < newdata_mmr60$cycle_end),"time_min"]), # purple
  				             col=rgb(red = 128, green = 0, blue = 128, alpha = 230, maxColorValue = 255),
  				             lty = "solid", lwd=4,lend="round",
  				             x1=newdata_mmr60$cycle_start, x2=newdata_mmr60$cycle_end);
  				  ablineclip(lm(d[which(d$time_min > newdata_mmr90$cycle_start &
  				                          d$time_min < newdata_mmr90$cycle_end),r]~
  				                  d[which(d$time_min > newdata_mmr90$cycle_start &
  				                            d$time_min < newdata_mmr90$cycle_end),"time_min"]), # orange
  				             col=rgb(red = 255, green = 127, blue = 80, alpha = 230, maxColorValue = 255),
  				             lty = "solid", lwd=4, lend="round",
  				             x1=newdata_mmr90$cycle_start, x2=newdata_mmr90$cycle_end);
  				  ablineclip(lm(d[which(d$time_min > newdata_mmr120$cycle_start &
  				                          d$time_min < newdata_mmr120$cycle_end),r]~
  				                  d[which(d$time_min > newdata_mmr120$cycle_start &
  				                            d$time_min < newdata_mmr120$cycle_end),"time_min"]),
  				             col=rgb(red = 0, green = 128, blue = 128, alpha = 230, maxColorValue = 255), # teal
  				             lty = "solid", lwd=4, lend="round",
  				             x1=newdata_mmr120$cycle_start, x2=newdata_mmr120$cycle_end);
				  }
				  # points(d[,r]~d$time_min, col="grey");
				  mtext(bquote(y == .(lm_coef[2])*x + .(lm_coef[1])), adj=1, padj=0, cex=0.8, line=0); # display equation
				  mtext(bquote(italic(R)^2 == .(format(r2, digits = 3))),adj=1, padj=0, cex=0.8, line=1)
				}
			}

			if (Ch_list[j]==4 & cycles>=4){

				DateTime_start<-as.character(d$DateTime[1])
				cycle_type<-"MMR"
				Ch<-paste("Ch",j, sep="")
				values_mmr<-as.data.frame(t(c(cycle_type, mmr_start, mmr_end, cycle_mmr, r2, m, b,
				                              round(min(d[,r_temp]), 4),
				                              round(max(d[,r_temp]), 4),
				                              round(mean(d[,r_temp]), 4),
				                              Ch, DateTime_start)))
				colnames(values_mmr)<-c("cycle_type","cycle_start","cycle_end",  "cycle_mmr", "r2" ,"m", "b" , "t_min", "t_max", "t_mean", "Ch", "DateTime_start")
				newdata_mmr<-rbind(newdata_mmr, values_mmr)

				newdata_mmr<-slidingSlope(d, Ch, data.MMR,r, r_temp, newdata_mmr, path)
				newdata_mmr60<-newdata_mmr[which(newdata_mmr$Ch == Ch & newdata_mmr$cycle_mmr == 60),]
				newdata_mmr90<-newdata_mmr[which(newdata_mmr$Ch == Ch & newdata_mmr$cycle_mmr == 90),]
				newdata_mmr120<-newdata_mmr[which(newdata_mmr$Ch == Ch & newdata_mmr$cycle_mmr == 120),]


				p4 %<a-% {
				  plot(d[,r]~time_min, data=d, col="grey", ylab=expression(paste(O2~(mg~L^{-1}))),xlab="Time (min)",main="MMR");
				  abline(lm(d[,r]~d$time_min), col="black",lwd=2);
				  if(nrow(newdata_mmr60)>0){
				    ablineclip(lm(d[which(d$time_min > newdata_mmr60$cycle_start &
  				                          d$time_min < newdata_mmr60$cycle_end),r]~
  				                  d[which(d$time_min > newdata_mmr60$cycle_start &
  				                            d$time_min < newdata_mmr60$cycle_end),"time_min"]), # purple
  				             col=rgb(red = 128, green = 0, blue = 128, alpha = 230, maxColorValue = 255),
  				             lty = "solid" , lwd=4,lend="round",
  				             x1=newdata_mmr60$cycle_start, x2=newdata_mmr60$cycle_end);
  				  ablineclip(lm(d[which(d$time_min > newdata_mmr90$cycle_start &
  				                          d$time_min < newdata_mmr90$cycle_end),r]~
  				                  d[which(d$time_min > newdata_mmr90$cycle_start &
  				                            d$time_min < newdata_mmr90$cycle_end),"time_min"]), # orange
  				             col=rgb(red = 255, green = 127, blue = 80, alpha = 230, maxColorValue = 255),
  				             lty = "solid", lwd=4, lend="round",
  				             x1=newdata_mmr90$cycle_start, x2=newdata_mmr90$cycle_end);
  				  ablineclip(lm(d[which(d$time_min > newdata_mmr120$cycle_start &
  				                          d$time_min < newdata_mmr120$cycle_end),r]~
  				                  d[which(d$time_min > newdata_mmr120$cycle_start &
  				                            d$time_min < newdata_mmr120$cycle_end),"time_min"]),
  				             col=rgb(red = 0, green = 128, blue = 128, alpha = 230, maxColorValue = 255), # teal
  				             lty = "solid", lwd=4, lend="round",
  				             x1=newdata_mmr120$cycle_start, x2=newdata_mmr120$cycle_end);
				  }
				  # points(d[,r]~d$time_min, col="grey");
				  mtext(bquote(y == .(lm_coef[2])*x + .(lm_coef[1])), adj=1, padj=0, cex=0.8, line=0); # display equation
				  mtext(bquote(italic(R)^2 == .(format(r2, digits = 3))),adj=1, padj=0, cex=0.8, line=1)
				}

			}


			#
			## non MMR plots
			if ((Ch_list[j]==1 & cycles==2) | (Ch_list[j]==1 & cycles==3) | (Ch_list[j]==1 & cycles>=4)){
			# if ((Ch_list[j]==1 & cycles>=2)){

				p2 %<a-% {
					plot(d2[,r]~time_min, data=d2, col="grey", ylab=expression(paste(O2~(mg~L^{-1}))),xlab="Time (min)",main="cycle2");
					abline(lm(d2[,r]~d2$time_min), col="red",lwd=2);
					mtext(bquote(y == .(lm_coef.2[2])*x + .(lm_coef[1])), adj=1, padj=0, cex=0.8, line=0); # display equation
					mtext(bquote(italic(R)^2 == .(format(r2.2, digits = 3))),adj=1, padj=0, cex=0.8, line=1);
				  if(!is.na(clean_cycle2)){
				    # print("here here")

				    if(start_c2==0 & end_c2 ==0){
				      legend("center","center", c("EXCLUDING \n from dataset"), col=c("red"), pch=" ", cex=2)
				    }else{
				      legend("topright", c("Time-adjusted"), col=c("red"), pch=" ", cex=1)
				    }
				  }
				}

				DateTime_start<-as.character(d2$DateTime[1])
				cycle_type<-"cycle2"
				Ch<-paste("Ch",j, sep="")

				if(is.na(clean_cycle2)){
  			  values_mmr<-as.data.frame(t(c(cycle_type, cycle_start[2], cycle_end[2],  cycle_mmr, r2.2, m.2, b.2,
  			                                round(min(d2[,r_temp]),4),
  			                                round(max(d2[,r_temp]), 4),
  			                                round(mean(d2[,r_temp]),4),
  			                                Ch, DateTime_start)))
  			  colnames(values_mmr)<-c("cycle_type","cycle_start","cycle_end",  "cycle_mmr", "r2" ,"m", "b" , "t_min", "t_max", "t_mean", "Ch", "DateTime_start")
  			  newdata_mmr<-rbind(newdata_mmr, values_mmr)
				}else{
				  if(start_c2==0 & end_c2 ==0){
				    message(paste("Take out cycle 2 from MMR trends - Channel: ", j))
				  }else{
            values_mmr<-as.data.frame(t(c(cycle_type, start_c2, end_c2,  cycle_mmr, r2.2, m.2, b.2,
                                          round(min(d2[,r_temp]),4),
                                          round(max(d2[,r_temp]), 4),
                                          round(mean(d2[,r_temp]),4),
                                          Ch, DateTime_start)))
            colnames(values_mmr)<-c("cycle_type","cycle_start","cycle_end",  "cycle_mmr", "r2" ,"m", "b" , "t_min", "t_max", "t_mean", "Ch", "DateTime_start")
            newdata_mmr<-rbind(newdata_mmr, values_mmr)
				  }
				}
			}

			# if ((Ch_list[j]==1 | Ch_list[j]==2 | Ch_list[j]==4) & cycles>=3){
			if ((Ch_list[j]==1 & cycles==3) | (Ch_list[j]== 1 & cycles==4) | (Ch_list[j]==2 & cycles==3) | (Ch_list[j]== 2 & cycles==4)){
				p3 %<a-% {
					plot(d3[,r]~time_min, data=d3, col="grey", ylab=expression(paste(O2~(mg~L^{-1}))),xlab="Time (min)",main="cycle3");
					abline(lm(d3[,r]~d3$time_min), col="red",lwd=2);
					mtext(bquote(y == .(lm_coef.3[2])*x + .(lm_coef.3[1])), adj=1, padj=0, cex=0.8, line=0); # display equation
					mtext(bquote(italic(R)^2 == .(format(r2.3, digits = 3))),adj=1, padj=0, cex=0.8, line=1);
				  if(!is.na(clean_cycle3)){

				    if(start_c3==0 & end_c3 ==0){
				      legend("center","center", c("EXCLUDING \n from dataset"), col=c("red"), pch=" ", cex=2)
				    }else{
				      legend("topright", c("Time-adjusted"), col=c("red"), pch=" ", cex=1)
            }
				  }
				}
				# print("here -here")
				# print(c(start_c3, end_c3))
				DateTime_start<-as.character(d3$DateTime[1])
				cycle_type<-"cycle3"
				Ch<-paste("Ch",j, sep="")
				if(is.na(clean_cycle3)){
          values_mmr<-as.data.frame(t(c(cycle_type, cycle_start[3], cycle_end[3],  cycle_mmr, r2.3, m.3, b.3,
                                        round(min(d3[,r_temp]),4),
                                        round(max(d3[,r_temp]), 4),
                                        round(mean(d3[,r_temp]),4),
                                        Ch, DateTime_start)))
          colnames(values_mmr)<-c("cycle_type","cycle_start","cycle_end",  "cycle_mmr", "r2" ,"m", "b" , "t_min", "t_max", "t_mean", "Ch", "DateTime_start")
          newdata_mmr<-rbind(newdata_mmr, values_mmr)
				}else{
				  if(start_c3==0 & end_c3 ==0){
				    message(paste("Take out cycle 3 from MMR trends - Channel: ", j))
				  }else{
				    values_mmr<-as.data.frame(t(c(cycle_type, start_c3, end_c3, cycle_mmr, r2.3, m.3, b.3,
				                                  round(min(d3[,r_temp]),4),
				                                  round(max(d3[,r_temp]), 4),
				                                  round(mean(d3[,r_temp]),4),
				                                  Ch, DateTime_start)))
				    colnames(values_mmr)<-c("cycle_type","cycle_start","cycle_end",  "cycle_mmr", "r2" ,"m", "b" , "t_min", "t_max", "t_mean", "Ch", "DateTime_start")
				    newdata_mmr<-rbind(newdata_mmr, values_mmr)
				  }
				}
			}

			# if ((Ch_list[j]==1 | Ch_list[j]==2 | Ch_list[j]==3) & cycles>=4){
			if ((Ch_list[j]==1 & cycles==4) | (Ch_list[j]== 2 & cycles==4) | (Ch_list[j]==3 & cycles==4)){

				p4 %<a-% {
					plot(d4[,r]~time_min, data=d4, col="grey", ylab=expression(paste(O2~(mg~L^{-1}))),xlab="Time (min)",main="cycle4");
					abline(lm(d4[,r]~d4$time_min), col="red",lwd=2);
					mtext(bquote(y == .(lm_coef.4[2])*x + .(lm_coef.4[1])), adj=1, padj=0, cex=0.8, line=0); # display equation
					mtext(bquote(italic(R)^2 == .(format(r2.4, digits = 3))),adj=1, padj=0, cex=0.8, line=1);
				  if(!is.na(clean_cycle4)){
				    if(start_c4==0 & end_c4 ==0){
				      legend("center","center", c("EXCLUDING \n from dataset"), col=c("red"), pch=" ", cex=2)
				    }else{
				      legend("topright", c("Time-adjusted"), col=c("red"), pch=" ", cex=1)
				    }
				  }
				}

				DateTime_start<-as.character(d4$DateTime[1])
				cycle_type<-"cycle4"
				Ch<-paste("Ch",j, sep="")

				if(is.na(clean_cycle4)){
				  values_mmr<-as.data.frame(t(c(cycle_type, cycle_start[4], cycle_end[4],  cycle_mmr, r2.4, m.4, b.4,
				                                round(min(d4[,r_temp]),4),
				                                round(max(d4[,r_temp]), 4),
				                                round(mean(d4[,r_temp]),4),
				                                Ch, DateTime_start)))
				  colnames(values_mmr)<-c("cycle_type","cycle_start","cycle_end",  "cycle_mmr", "r2" ,"m", "b" , "t_min", "t_max", "t_mean", "Ch", "DateTime_start")
				  newdata_mmr<-rbind(newdata_mmr, values_mmr)
				}else{
				  if(start_c4==0 & end_c4 ==0){
				    message(paste("Take out cycle 4 from MMR trends - Channel: ", j))
				  }else{
				    values_mmr<-as.data.frame(t(c(cycle_type, start_c4, end_c4,  cycle_mmr, r2.4, m.4, b.4,
				                                  round(min(d4[,r_temp]),4),
				                                  round(max(d4[,r_temp]), 4),
				                                  round(mean(d4[,r_temp]),4),
				                                  Ch, DateTime_start)))
				    colnames(values_mmr)<-c("cycle_type","cycle_start","cycle_end",  "cycle_mmr", "r2" ,"m", "b" , "t_min", "t_max", "t_mean", "Ch", "DateTime_start")
				    newdata_mmr<-rbind(newdata_mmr, values_mmr)
				  }
				}
			}


			# null plot
			p_null %<a-% {
				plot(1, type="n", axes=T, xlab="", ylab="");
				text(1, 1, "MO2 NOT MEASURED",cex = .8)
			}

			png(plotname, width=20, height=5, units="in",  res=200)
			par(mfrow=c(1,5))
				if (cycles==1){
					p
					p1
					p_null
					p_null
					p_null
				}

				if (cycles==2){
					if (Ch_list[j]==1){
						p
						p1
						p2
						p_null
						p_null
					}
					if (Ch_list[j]==2){
						p
						p_null
						p2
						p_null
						p_null
					}
				}

				if (cycles==3){
					if (Ch_list[j]==1){
						p
						p1 # MMR
						p2
						p3
						p_null
					}
					if (Ch_list[j]==2){
						p
						p_null
						p2 # MMR
						p3
						p_null
					}
					if (Ch_list[j]==3){
						p
						p_null
						p_null
						p3 # MMR
						p_null
					}
				}

				if (cycles>=4){
					if (Ch_list[j]==1){
						p
						p1 # MMR
						p2
						p3
						p4
					}
					if (Ch_list[j]==2){
						p
						p_null
						p2 # MMR
						p3
						p4
					}
					if (Ch_list[j]==3){
						p
						p_null
						p_null
						p3 # MMR
						p4
					}
					if (Ch_list[j]==4){
						p
						p_null
						p_null
						p_null
						p4 # MMR
					}
				}

			dev.off()
			#
			#  save file now -- write code
		}else{
				message(paste("NO Channel - ",j, sep=""))
		}

	#~ 	assign("newdata_mmr", newdata_mmr, envir=.GlobalEnv)

		# Channel_mmr<-newdata_mmr
		return(newdata_mmr)

	}# end of Channel_mmr function



	newdata_mmr<-matrix(ncol=12,nrow=0)
	colnames(newdata_mmr)<-c("cycle_type", "cycle_start","cycle_end",  "cycle_mmr", "r2" ,"m", "b" , "t_min", "t_max", "t_mean", "Ch", "DateTime_start")


	if (path == "."){
    filename<-paste( gsub('.{4}$', '', data.MMR), "_analyzed.csv", sep='')
    filename2<-paste( gsub('.{4}$', '', data.MMR), "_analyzed.csv", sep='') # save in the EPOC_AS etc folder for final MMR and full EPOC analysis

	}else{
	  filename<-paste("../csv_analyzed/", gsub('.{4}$', '', data.MMR), "_analyzed.csv", sep='')
    filename2<-paste("../../MMR_SMR_AS_EPOC/csv_input_files/", gsub('.{4}$', '', data.MMR), "_analyzed.csv", sep='') # save in the EPOC_AS etc folder for final MMR and full EPOC analysis

	}


	dataMMR<-read.csv(data.MMR)
	if(as.character(dataMMR$Ch1_O2[1])=="--- " || as.character(dataMMR$Ch1_O2[1])=="---"){
		dataMMR$Ch1_O2<-0
	}
	if(as.character(dataMMR$Ch2_O2[1])=="--- " || as.character(dataMMR$Ch2_O2[1])=="---"){
		dataMMR$Ch2_O2<-0
	}
	if(as.character(dataMMR$Ch3_O2[1])=="--- " || as.character(dataMMR$Ch3_O2[1])=="---"){
		dataMMR$Ch3_O2<-0
	}
	if(as.character(dataMMR$Ch4_O2[1])=="--- " || as.character(dataMMR$Ch4_O2[1])=="---"){
		dataMMR$Ch4_O2<-0
	}

	n<-length(which(!dataMMR[1,c(4, 6:8)]==0))


	dataMMR$date<-as.character(dataMMR$date)
	dataMMR$time<-as.character(dataMMR$time)

	if (date_format== "m/d/y"){
  	DateTime<- chron(dates.=dataMMR$date,times.=dataMMR$time,format=c('m/d/y','h:m:s')) # Eliason Lab firesting date format
#~ 	DateTime<- chron(dates=dataMMR$date,times=dataMMR$time,format=c('y-m-d','h:m:s')) # Healy firesting date format
	}
	if (date_format== "d/m/y"){
    # 8 channel firesting with date (day first) European style
    DateTime<- chron(dates.=dataMMR$date,times.=dataMMR$time,format=c('d/m/y','h:m:s')) # Eliason Lab firesting date format
	}
	if (date_format== "y-m-d"){
    # 8 channel firesting with date (day first) European style
    DateTime<- chron(dates.=dataMMR$date,times.=dataMMR$time,format=c('Y-m-d','h:m:s')) # Eliason Lab firesting date format
  }

	dataMMR$DateTime<-DateTime

	dataMMR$hr<-round(dataMMR$time_sec/3600,2)
	dataMMR$hr2<-round(dataMMR$time_sec/3600,0)
	dataMMR$time_min<-round(dataMMR$time_sec/60,2)

	dataMMR$Ch1_temp<-as.numeric(as.character(dataMMR$Ch1_temp))

	if (N_Ch==8){
  	dataMMR$Ch2_temp<-as.numeric(as.character(dataMMR$Ch2_temp))
  	dataMMR$Ch3_temp<-as.numeric(as.character(dataMMR$Ch3_temp))
  	dataMMR$Ch4_temp<-as.numeric(as.character(dataMMR$Ch4_temp))
  		#
    temp_mean2<-mean(dataMMR$Ch2_temp, na.rm=TRUE)
    temp_max2<-max(dataMMR$Ch2_temp, na.rm=TRUE)
    temp_min2<-min(dataMMR$Ch2_temp, na.rm=TRUE)

    temp_mean3<-mean(dataMMR$Ch3_temp, na.rm=TRUE)
    temp_max3<-max(dataMMR$Ch3_temp, na.rm=TRUE)
    temp_min3<-min(dataMMR$Ch3_temp, na.rm=TRUE)

    temp_mean4<-mean(dataMMR$Ch4_temp, na.rm=TRUE)
    temp_max4<-max(dataMMR$Ch4_temp, na.rm=TRUE)
    temp_min4<-min(dataMMR$Ch4_temp, na.rm=TRUE)
  }

	dataMMR$Ch1_O2<-as.numeric(as.character(dataMMR$Ch1_O2))
	dataMMR$Ch2_O2<-as.numeric(as.character(dataMMR$Ch2_O2))
	dataMMR$Ch3_O2<-as.numeric(as.character(dataMMR$Ch3_O2))
	dataMMR$Ch4_O2<-as.numeric(as.character(dataMMR$Ch4_O2))

	#temperature
	temp_mean1<-mean(dataMMR$Ch1_temp, na.rm=TRUE)
	temp_max1<-max(dataMMR$Ch1_temp, na.rm=TRUE)
	temp_min1<-min(dataMMR$Ch1_temp, na.rm=TRUE)

	Ch_list<-c(mmr_Ch1, mmr_Ch2, mmr_Ch3, mmr_Ch4) ## specifies in which cycle does the channel has MMR at
	rows<-c(4,6,7,8) # the order Ch1, Ch2, Ch3, Ch4

	if(N_Ch==8){
  	rows_temp<-c(5,9,10,11) # the order Ch1, Ch2, Ch3, Ch4
	  }else {
	 	rows_temp<-c(5,5,5,5) # the order Ch1, Ch2, Ch3, Ch4
	}

	if(clean_Ch1[1] == 0 & mmr_Ch1!=0){
    clean_Ch1[1] <- cycle_start[mmr_Ch1]
	}
	if(clean_Ch2[1] == 0 & mmr_Ch2!=0){
    clean_Ch2[1] <- cycle_start[mmr_Ch2]
	}
	if(clean_Ch3[1] == 0 & mmr_Ch3!=0){
    clean_Ch3[1] <- cycle_start[mmr_Ch3]
	}
	if(clean_Ch4[1] == 0 & mmr_Ch4!=0){
    clean_Ch4[1] <- cycle_start[mmr_Ch4]
	}

	if(clean_Ch1[2] == 0 & mmr_Ch1!=0){
    clean_Ch1[2] <- cycle_end[mmr_Ch1]
	}
	if(clean_Ch2[2] == 0 & mmr_Ch2!=0){
    clean_Ch2[2] <- cycle_end[mmr_Ch2]
	}
	if(clean_Ch3[2] == 0 & mmr_Ch3!=0){
    clean_Ch3[2] <- cycle_end[mmr_Ch3]
	}
	if(clean_Ch4[2] == 0 & mmr_Ch4!=0){
    clean_Ch4[2] <- cycle_end[mmr_Ch4]
	}

  clean_start_mmr<-c(clean_Ch1[1], clean_Ch2[1], clean_Ch3[1], clean_Ch4[1]) # cleaning start
  # cycle_start
  clean_end_mmr<-c(clean_Ch1[2], clean_Ch2[2], clean_Ch3[2], clean_Ch4[2]) # cleaning end
  # cycle_end

  # if the last cycle end is provided to be longer (min) than the actual file is, then rewrite the end time to be equal to the time the respo is run.
	if (cycle_end[length(cycle_end)]>dataMMR$time_min[nrow(dataMMR)]){
		cycle_end[length(cycle_end)]<-dataMMR$time_min[nrow(dataMMR)]
	}



      Box_n_data<-(which(!is.na(str_locate(data.MMR, c("box", "BOX", "Box"))))[2])

      if(!is.na(inv.data)){
        data2<-read.csv(inv.data)

        # print(data2)
      }

      	for (j in 1:4){


      	  if(!is.na(inv.data)){


        	   if(j==1){inv.data.clean<-data2[which(grepl(substr(data.MMR, start=1, stop=5), as.character(data2$date)) & as.numeric(data2$channel)==1 &
        			                        grepl(substr(data.MMR, start = (str_locate(data.MMR, c("box", "BOX", "Box"))[Box_n_data])+1,
        			                                     stop = (str_locate(data.MMR, c("box", "BOX", "Box"))[Box_n_data])+1), as.character(data2$box))),]
        	   }
        	   if(j==2){inv.data.clean<-data2[which(grepl(substr(data.MMR, start=1, stop=5), as.character(data2$date)) & as.numeric(data2$channel)==2 &
        			                        grepl(substr(data.MMR, start = (str_locate(data.MMR, c("box", "BOX", "Box"))[Box_n_data])+1,
        			                                     stop = (str_locate(data.MMR, c("box", "BOX", "Box"))[Box_n_data])+1), as.character(data2$box))),]
        	   }
        	   if(j==3){inv.data.clean<-data2[which(grepl(substr(data.MMR, start=1, stop=5), as.character(data2$date)) & as.numeric(data2$channel)==3 &
        			                        grepl(substr(data.MMR, start = (str_locate(data.MMR, c("box", "BOX", "Box"))[Box_n_data])+1,
        			                                     stop = (str_locate(data.MMR, c("box", "BOX", "Box"))[Box_n_data])+1), as.character(data2$box))),]
        	   }
        	   if(j==4){inv.data.clean<-data2[which(grepl(substr(data.MMR, start=1, stop=5), as.character(data2$date)) & as.numeric(data2$channel)==4 &
        			                        grepl(substr(data.MMR, start = (str_locate(data.MMR, c("box", "BOX", "Box"))[Box_n_data])+1,
        			                                     stop = (str_locate(data.MMR, c("box", "BOX", "Box"))[Box_n_data])+1), as.character(data2$box))),]
        	   }

      	  }else{

      	    inv.data.clean<-as.data.frame(matrix(nrow=0, ncol=0))
      	  }

      		newdata_mmr<-Channel_mmr(data.MMR, dataMMR, cycle_start, cycle_end, clean_start_mmr, clean_end_mmr, Ch_list, cycles, j, rows, rows_temp, newdata_mmr, path, inv.data.clean)
      	}


	newdata_mmr$m<-abs(as.numeric(as.character(newdata_mmr$m)))

	write.csv(file=filename, newdata_mmr, row.names=FALSE)
	write.csv(file=filename2, newdata_mmr, row.names=FALSE)

#~ 	assign("newdata_mmr", newdata_mmr, envir=.GlobalEnv)

}

