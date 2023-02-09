#' Title
#'
#' @param AnimalID Indicates individual ID; must be a vector of 4 characters. When missing, enter "NA"
#' @param BW.animal Indicates individual mass; must be a vector of 4 characters. When missing, enter "0"
#' @param resp.V Indicated the volume (L) of respirometer chambers; must be a vector of 4 numbers (e.g., c(1, 1, 1, 1), for four 1-L respirometers)
#' @param r2_threshold_smr R2 threshold for SMR, measurements below the threshold are excluded
#' @param sda_threshold_level Indicates the threshold relevant to an individual’s SMR to calculate the end time of recovery (the time at which metabolic rate has returned to sda_threshold). The default is the SMR level (1; At 100 percent SMR). To use 120 percent SMR as threshold for full digestion, enter sda_threshold = 1.2.
#' @param scaling_exponent_smr Body mass scaling exponent to correct SMR values for body size. MR=aBM^b (MR = metabolic rate, BW = body mass, a = scaling coefficient [the intercept], and b = scaling exponent [the power term])
#' @param date_format The date format used in the original FireSting data files. Must specify one of the following: "m/d/y", "d/m/y", "y-m-d"
#' @param data.SDA The name of the SDA data file (“…analyzed.csv”; a character string); an output file from the SMR function.
#' @param analyzed_MR The name of the data file with estimated SMR values
#' @param SMR_calc Logical, indicates whether to calculate SMR from the present data trend
#' @param SMR_vals allows a user to input their own SMR values
#' @param N_Ch The number of channels of the FireSting. It must be either 4 or 8.
#' @param drop_ch Indicates which channel is dropped or entirely excluded from the analysis. Must be a numerical vector, e.g., c(1,3)
#' @param end_SDA_Ch Manually assigned end time (min) of digestion (SDA). The time can be assigned for each channel independently; use a numerical vector of 4 variables, one for each channel (e.g., c(120, 120, NA, 180), for 2 h, 2h, SMR level, and 3 h SDA end times, respectively)
#' @param MLND Logical argument. If TRUE, SMR is estimated also using Mean Lowest Normal Distribution analysis. More details in Chabot et al 2016.
#' @param verbose.MLND From MLND: A logical controlling if a text progress bar from MLND is displayed during the fitting procedure. (see 'verbose' in mclust package functions).
#' @param background_prior The name of the analyzed background “…analyzed.csv” file, an output file from the SMR function (respiration in an empty respirometer measured before the respirometry trial).
#' @param background_post The same as background_prior, only post respirometry trial
#' @param background_slope Manually assigned background slope used to correct metabolic rate in all individuals. Provide numeric value in mgO2 L-1 h-1v
#' @param background.V Manually assigned respirometer volumes (L). A vector with 4 numeric variables, one for each channel.
#' @param background_gr Specify whether to assume that bacterial growth (thus respiration rates) changed linearly or in exponentially across the duration of the respirometry trial. Must specify either "linear" or "exp": metabolic rate values across the given trial are corrected using the estimated background values from the indicated growth curve. Both background_prior and background_post must be provided to enable this.
#' @param match_background_Ch Logical. If TRUE, the background respiration is estimated and applied channel-specific. The background_prior and background_post are used to estimate background respiration specific to each channel (or individual), which then is used to correct each individual’s MO2 independently. By default, the mean background respiration rate from all channels is calculated and applied to correct all individual’s MO2
#' @param local_path Logical. If TRUE (default) all returned files will be saved in the local working directory.
#' @param handling_delay a delay of SDA due to handling, anesthsia, or any other experimental factor. in hours
#' @param begin_smr_hr_zero logical. Wether to force the integral SDA costs still be calculated from hour 0, that is assumed to be at SMR level (refer to sda_threshold_level)
#'
#' @importFrom stats lm coef var integrate predict quantile sd smooth.spline complete.cases
#' @import graphics
#' @import grDevices
#' @importFrom gridExtra grid.arrange
#' @import scales
#' @import ggplot2
#' @import utils
#' @importFrom dplyr filter top_n arrange
#' @importFrom tidyr spread
#' @importFrom mclust Mclust densityMclust MclustBootstrap
#' @importFrom DescTools AUC
#' @importFrom magrittr %>%
#'
#' @return The output from \code{\link{print}}
#' @export
#'
SDA<-function(AnimalID,
              BW.animal,
              resp.V,
              r2_threshold_smr,
              sda_threshold_level,
              scaling_exponent_smr = 1,
              date_format = c("%Y-%m-%d %H:%M:%S", "%m/%d/%Y %H:%M:%S"),
              data.SDA = NULL,
              analyzed_MR = NULL,
              SMR_calc = TRUE,
              SMR_vals = c(NULL, NULL, NULL, NULL),
              drop_ch = NULL,
              N_Ch = 4,
              end_SDA_Ch = NA,
              MLND = TRUE,
              verbose.MLND = FALSE,
              background_prior = NULL,
              background_post = NULL,
              background_slope = NULL,
              background.V = NULL,
              match_background_Ch = FALSE,
              background_gr = NULL,
              local_path = TRUE,
              handling_delay = c(0,0,0,0),
              begin_smr_hr_zero=FALSE){

  # *************** SETUP START
  DateTime_start <- back_regression1 <- back_regression2 <- back_regression3 <- back_regression4 <- NULL
  bw <- hour <- m <- min_start <- mo2_mean <- mo2_min <- mo2_perc <- quantiles <- smr_method <- smr_val <- NULL


  if(!length(as.vector(date_format))==2){
    stop_function<-TRUE
    if(stop_function) {
      stop("Cannot interpret 'date_format' format. \n Please provide a vector of two stating: i) the date time format and ii) timezone. \n Default is: c(\"%Y-%m-%d %H:%M:%S\", \"GMT\"). Argument is passed to strptime() ")
    }
  }

  # read in files SDA ------
  if(file.exists(data.SDA[1]) | file.exists(paste("./SDA/csv_input_files/", data.SDA[1], sep=""))){ # after running through RMRrepeat - this will be saved in csv input files

    if(file.exists(paste("./SDA/csv_input_files/", data.SDA[1], sep=""))){
      data_glued<-read.csv(paste("./SDA/csv_input_files/", data.SDA[1], sep=""))
      if(length(data.SDA) > 1){
        message("Multiple SDA files are combined")
        for(i in 2:length(data.SDA)){
          data_glued0 <- read.csv(paste("./SDA/csv_input_files/", data.SDA[i], sep=""))
          data_glued <- rbind(data_glued, data_glued0)
        }
      }
    }

    if(file.exists(data.SDA[1])){
      data_glued<-read.csv(data.SDA[1])
      if(length(data.SDA) > 1){
        message("Multiple SDA files are combined")
        for(i in 2:length(data.SDA)){
          data_glued0 <- read.csv(data.SDA[i])
          data_glued <- rbind(data_glued, data_glued0)
        }
      }
    }

  }else{
    stop_function<-TRUE
    if(stop_function){
      stop("Cannot locate the indicated data.SDA data file.")
    }
  }

  # first order all data by channels
  data_glued <- data_glued[order(data_glued$Ch), ] # even if only one file

  if(!is.na(strptime(data_glued$DateTime_start[1], format = date_format[1], tz = ""))){
    data_glued$DateTime_start<- strptime(data_glued$DateTime_start, format = date_format[1], tz = "")
  }
  if(!is.na(strptime(data_glued$DateTime_start[1], format = date_format[2], tz = ""))){
    data_glued$DateTime_start<- strptime(data_glued$DateTime_start, format = date_format[2], tz = "")
  }

  if(length(data.SDA) > 1){
    data_glued$time_diff<-NA

    for(i in 2:nrow(data_glued)){
      data_glued$time_diff[i]<-data_glued$min_start[1]+(as.numeric(difftime(data_glued$DateTime_start[i],data_glued$DateTime_start[1], units="min")))
    }

    data_glued$time_diff[1]<-data_glued$min_start[1]
  	data_glued$min_start<-data_glued$time_diff
  	data_glued<-data_glued[,1:16]
  }
	d_SMR<-data_glued


# end read in file
  SDA.spar<-function(spar,d, SDAdata, b, sda_threshold, end_SDA, begin_smr_hr_zero){

    sda_threshold<-as.numeric(sda_threshold)
    d$hour<-as.numeric(as.character(d$hour))

    d<-d[complete.cases(d),]
    fit<-smooth.spline(d$hour,d$mo2_min, spar=spar)
    f = function(x) {predict(fit, x)$y}

    end<-round(d$hour[nrow(d)],1)
    newx<-seq(0,end, by=0.1)
    newy<-predict(fit, newx, deriv=0)
    newy<-lapply(newy, as.numeric)
    newVal<-data.frame(Reduce(cbind,newy)) # creating dummy dataset with predicted values

    SDAdata.temp<-matrix(ncol=15, nrow=0)
    colnames(SDAdata)<-c("ID","SMR","spar", "SDA_integrated", "end_SDA_estim_hr", "SMR_intergrated", "peak_SDA", "time_peak_SDA", "percentSMR_peak_SDA", "MO2_SDA_full", "peak_SDA_max", "time_peak_SDA_max", "peak_SDA_mean", "time_peak_SDA_mean", "smr_type")

    #2-1- establish SMR threshold m1 and b1 for SMR (b1= y intercept, m1=slope)
    m0 <- 0 # slope
    f.smr <- function(x)(m0*x)+b # smr function

  # Calculate the area under the SMR for all types of SMR calculated in the MMr_SMR_analyze function
  #2-2 the EPOC cuttoff in the smoothed function / this is used also for the SMR block
  # providing end_SDA value manually (in minutes)
  if(is.na(end_SDA)){
    if (begin_smr_hr_zero==TRUE){
      end_SDA <- newVal$init[(which(round(newVal$V2[2:nrow(newVal)], 3)<=b))[1]] # force the function to look beyond te first value, which is set to SMR level
    }else{
      end_SDA <- newVal$init[(which(round(newVal$V2, 3)<=b))[1]]
    }
    if(is.na(end_SDA)){
      end_SDA<-newVal$init[nrow(newVal)]
      message(paste("Smooth level spar:", spar, "MR does not reach the chosen SMR levels: ",b, ": end tie of SDA is the end of the trial"))
    }
  }


  # if(end_SDA > max(newVal$init)){
  #   message("The SDA end value is greater than the duration of respirometry trial | end_SDA set to the last value")
  #   end_SDA<-newVal$init[nrow(newVal)]
  # }



  #3 integrate recovery curve with to the provided end EPOC time
  f.int = function(x) {integrate(f, lower=0, upper=end_SDA)$value}
  f.int.vec = Vectorize(f.int, vectorize.args='x')
  f.int.vec = Vectorize(f.int, vectorize.args='x')
  full<-f.int.vec(end_SDA)
  # SMR block
  SMR<-integrate(f.smr, lower=0, upper=end_SDA)$value
  SDA_full<-round(full-SMR,3) # the costs of digestion
  MO2_full<-round(newVal$V2[which(newVal$init==end_SDA)],3)

  # end_SDA_row<-which(round(d$hour)==round(end_SDA))[1]
  end_SDA_row <- which(abs(d$hour - end_SDA) == min(abs(d$hour - end_SDA)))[1]
  peak_SDA<-max(d$mo2_min[1:end_SDA_row])[1]
  peak_SDA_max<-max(d$mo2_max[1:end_SDA_row])[1]
  peak_SDA_mean<-max(d$mo2_mean[1:end_SDA_row])[1]

  time_peak_SDA<-d$hour[which(d$mo2_min==peak_SDA)]
  time_peak_SDA_max<-d$hour[which(d$mo2_max==peak_SDA_max)[1]]
  time_peak_SDA_mean<-d$hour[which(d$mo2_mean==peak_SDA_mean)[1]] # in h not saved in the dataframe
  percentSMR_peak_SDA<-round((peak_SDA/b)*100,2)

  values<-as.data.frame(t(c(as.character(d$ID[1]),
                            b, spar, SDA_full, end_SDA,
                            SMR, peak_SDA, time_peak_SDA, percentSMR_peak_SDA, MO2_full,
                            peak_SDA_max, time_peak_SDA_max, peak_SDA_mean, time_peak_SDA_mean,  sda_threshold)))


  colnames(values)<-c("ID","SMR","spar", "SDA_integrated", "end_SDA_estim_hr",
                      "SMR_intergrated", "peak_SDA", "time_peak_SDA", "percentSMR_peak_SDA", "MO2_SDA_full",
                      "peak_SDA_max", "time_peak_SDA_max", "peak_SDA_mean", "time_peak_SDA_mean", "smr_type")

  SDAdata.temp<-rbind(SDAdata.temp,values)
  scale<-c(0.02,0.07,0.12,0.17,0.22)

  plot(x=range(newVal$init), xlim=c(0,max(d$hour, na.rm=TRUE)+2), ylim=c(0,max(d$mo2_min, na.rm=TRUE)+2),type='n', ylab="MO2", xlab="time (h)")
  points(d$hour,d$mo2_min, main=spar, col="grey50", pch=21)
  lines(newy[[1]], newy[[2]], col="black")
  abline(h=b, col="red",lty=1, lwd=1)
  abline(v=end_SDA, col="red", lty=2)
  legend("topright", legend=paste("spar = ", spar, "\n blue diam: peak SDA mean \n red diam: peak SDA" ))
  # text(x=(max(d$mo2_min)+50)-(0.01*(max(d$hour)+50)),y=(max(d$mo2_min)+2)-(scale[i]*(max(d$mo2_min)+2)), label=paste("EPOC=",EPOC_full,"/ ",smr_type, sep=""), cex=0.8, col=col_smr[i], pos=2)
  points(x=time_peak_SDA, y=peak_SDA, pch=23, col="red",cex=2)
  points(x=time_peak_SDA_mean, y=peak_SDA_mean, pch=23, col="blue",cex=2)
  # text(x=(max(d$time_mo2)+50)-(0.01*(max(d$time_mo2)+50)),y=(max(d$mo2)+2)-(scale[i]*(max(d$mo2)+2)), label=paste("EPOC=",EPOC_full,"/ ",smr_type, sep=""), cex=0.8, col=col_smr[i], pos=2)


  SDAdata<-rbind(SDAdata,SDAdata.temp)

  return(SDAdata)

} # end SDA.spar

  graphics.off()

  filename.SMR<-paste(gsub('.{4}$', '',data.SDA[1]), "_SMR", sep="")

  if (local_path | !dir.exists("SDA")){
		plotname.freq<-paste( filename.SMR,"_PLOT_SMR_analyses.png", sep="")
		plotname.smr.meth<-paste( filename.SMR,"_PLOT_SMR_methodsALL.png", sep="")
	}else{
		plotname.freq<-paste("./SDA/plots_min_values_SMR/", filename.SMR,"_PLOT_SMR_analyses.png", sep="")
		plotname.smr.meth<-paste("./SDA/plots_methods_sum_SMR/", filename.SMR,"_PLOT_SMR_methodsALL.png", sep="")
  }

  newdata.smr<-as.data.frame(matrix(ncol=17, nrow=0))
  names(newdata.smr)<-c("filename", "ID", "Ch", "BW","t_min","t_max", "t_mean", "N_mo2", #8
                        "smr_mean10minVal","smr_SD10minVal", "smr_CV10minVal", "SMR_low10quant","SMR_low15quant","SMR_low20quant", #6
                        "smr_mlnd", "smr_CVmlnd", "smr_Nmlnd")#3

  cols = c(4:17)
  newdata.smr[,cols] <- lapply(newdata.smr[,cols], as.character)
  newdata.smr[,cols] <- lapply(newdata.smr[,cols], as.numeric)
  cols2 = c(1:3)
  newdata.smr[,cols2] <- lapply(newdata.smr[,cols2], as.character)

  # *************** SETUP END





  # background  -----
  if((is.null(background.V) & !is.null(background_slope)) |(!is.null(background.V) & is.null(background_slope))) {
    stop_function <- TRUE
    if(stop_function) stop("Must provide background unique slope and volume together, or provide a background file with same volumes as animal respirometers")
  }

	if(c(is.null(background_prior) & is.null(background_post)) & match_background_Ch){
    stop_function <- TRUE
	  if(stop_function) stop("If match_background_Ch = TRUE, must provide at least one datafiles, either background measured before or after the respirometry trial")
	}

  # 1. find what channels recorded background
  if (!is.null(background_post) | !is.null(background_prior)) {

    if(!is.null(background_prior)){

      if(file.exists(background_prior) | file.exists(paste("./SDA/csv_input_files/", background_prior, sep=""))){ # after running through RMRrepeat - this will be saved in csv input files
      	  if(file.exists(paste("./SDA/csv_input_files/", background_prior, sep=""))){
            back_prior<-read.csv(paste("./SDA/csv_input_files/", background_prior, sep=""))
          }
          if(file.exists(background_prior)){
            back_prior<-read.csv(background_prior)

          }
      	}else{
          stop_function<-TRUE
          if(stop_function){
            stop("Cannot locate the indicated background_prior data file.")
          }
      }
    }

    if(!is.null(background_post)){
      if(file.exists(background_post) | file.exists(paste("./SDA/csv_input_files/", background_post, sep=""))){ # after running through RMRrepeat - this will be saved in csv input files
    	  if(file.exists(paste("./SDA/csv_input_files/", background_post, sep=""))){
          back_post<-read.csv(paste("./SDA/csv_input_files/", background_post, sep=""))
        }
        if(file.exists(background_post)){
          back_post<-read.csv(background_post)
        }
    	}else{
        stop_function<-TRUE
        if(stop_function){
          stop("Cannot locate the indicated background_post data file.")
        }
      }
    }


    # Jan 4 2020: make linear regression over time and correct background based on a predicted value
    # create
    if (!is.null(background_gr)){
      if(!exists("back_prior") | !exists("back_post")){
        stop("Missing a file to model the growth of bacteria, must provide both: \"back_prior\" \"and back_post\" files")
      }

      if (!dir.exists("plots_background") | !dir.exists("./SDA/plots_background")){
        if(local_path & !dir.exists("plots_background")){
          dir.create(file.path("./plots_background"), recursive = TRUE)
        }
        if(local_path == FALSE &!dir.exists("./SDA/plots_background")){
          dir.create(file.path("./SDA/plots_background"), recursive = TRUE)
        }
      }

    if(!is.na(strptime(back_prior$DateTime_start[1], format = date_format[1], tz = ""))){
      back_prior$DateTime_start<- strptime(back_prior$DateTime_start, format = date_format[1], tz = "")
  	  back_post$DateTime_start<- strptime(back_post$DateTime_start, format = date_format[1], tz = "")
    }
    if(!is.na(strptime(back_prior$DateTime_start[1], format = date_format[2], tz = ""))){
      back_prior$DateTime_start<- strptime(back_prior$DateTime_start, format = date_format[2], tz = "")
  	  back_post$DateTime_start<- strptime(back_post$DateTime_start, format = date_format[2], tz = "")
    }

      back_all<- rbind(back_prior, back_post)
      back_all$DateTime_start<- as.POSIXct(back_all$DateTime_start)

      # back_regression_plot<-ggplot(data=back_all, aes(x=DateTime_start, y=m, colour=Ch, group=Ch))+
      #   geom_point()+
      #   geom_smooth(method="lm", se=FALSE)+
      #   theme_bw()+
      #   theme(axis.text.x = element_text(angle = 45))+
      #   facet_grid(Ch~.)
      # png(plotname.backgr.analysis, width=4, height=7, res=300, units="in")
      #   print(back_regression_plot)
      # dev.off()


      if (match_background_Ch==TRUE){
        back_ch_regressions<-list()

        back_ch<-length(unique(back_all$Ch))

        for(i in 1:back_ch){
          back_ch_d<- back_all[back_all$Ch==(unique(back_all$Ch))[i],]
          Ch<-substr(as.character(back_ch_d$Ch[1]), start=3, stop=3)

          back_regression_name<- paste("back_regression", Ch, sep="") # channel names with a channel # at the end

          if(background_gr == "linear"){
            regression<- lm(m~DateTime_start, data = back_ch_d)
            if(i==1){
              message("Background: linear change of bacterial respiration. Regressions are specific to the channel (respirometer)")
            }
          }
          if(background_gr == "exp"){
            # print(head(back_ch_d))
            regression<- lm(log(m)~DateTime_start, data = back_ch_d)
            if(i==1){
              message("Background: exponential change of bacterial respiration. Regressions are specific to the channel (respirometer)")
            }
          }

          # regression<- lm(m~DateTime_start, data = back_ch_d)
          assign(back_regression_name, regression)
          back_ch_regressions[[i]] <- assign(back_regression_name, regression)

        }
        # WORK
        # save background regressions

      }else{# end for ch specific regressions
        # get one background slope based on all datapoints collected
        if(background_gr == "linear"){
          back_regression<- lm(m~DateTime_start, data = back_all)
          message("Background: linear change of bacterial respiration rates. One regression for all channels")
        }
        if(background_gr == "exp"){
          back_regression<- lm(log(m)~DateTime_start, data = back_all)
          message("Background: exponential change of bacterial respiration rates. One regression for all channels")
        }
      }


    } # end for getting linear regressions for the background

    if (is.null(background_gr)){ # no growth of bacteria
      if (!is.null(background_prior)){
        back_ch_prior<-list()
        back_ch_prior_names<-list()

        back_ch<-length(unique(back_prior$Ch))

        for( i in 1:back_ch){
          back_ch_d<-back_prior[back_prior$Ch==(unique(back_prior$Ch))[i],]
          Ch<-substr(as.character(back_ch_d$Ch[1]), start=3, stop=3)

          back_m_name<-paste("back_m_prior", Ch, sep="") # channel names with a channel # at the end
          mean_m<-sum(back_ch_d$m)/nrow(back_ch_d) # get average slopes from the prior background cycles can be limitless essentially
          assign(back_m_name, mean_m)
          back_ch_prior[[i]] <- assign(back_m_name, mean_m)
          back_ch_prior_names[[i]] <- back_m_name

        }

      }
      # 3. estiamte one background slope mean to be used in MR corrections
      if (!is.null(background_post)){

        back_ch_post<-list()
        back_ch_post_names<-list()
        back_ch<-length(unique(back_post$Ch))

        for( i in 1:back_ch){

          back_ch_d<-back_post[back_post$Ch==(unique(back_post$Ch))[i],]
          Ch<-substr(as.character(back_ch_d$Ch[1]), start=3, stop=3)

          back_m_name<-paste("back_m_post", Ch, sep="")
          mean_m<-sum(back_ch_d$m)/nrow(back_ch_d)

          assign(back_m_name, mean_m)

          back_ch_post[[i]] <- assign(back_m_name, mean_m)
          back_ch_post_names[[i]] <- back_m_name

        }
      }

      # get a list of our values
      # should look for possible variables:
      # back_m_post1
      # back_m_post2
      # back_m_post3
      # back_m_post4
      # back_m_prior1
      # back_m_prior2
      # back_m_prior3
      # back_m_prior4

      if (match_background_Ch==TRUE){
        # conditions possible:
        # prior only
        # post only
        # prior and post for a specific channel

        if (!is.na(background_prior)){
          ch_available<-as.numeric(substr(as.character(unique(back_prior$Ch)),start=3, stop=3))
        }
        if (!is.na(background_post)){
          ch_available<-as.numeric(substr(as.character(unique(back_post$Ch)),start=3, stop=3))
        }

        for (i in 1:back_ch){

          # prior only condition
          if(exists(paste("back_m_prior", ch_available[i], sep="")) & !exists(paste("back_m_post", ch_available[i], sep=""))){

            if(back_ch_prior_names[[i]]==paste("back_m_prior", ch_available[i], sep="")){
              message("matching background Channels")

              back_m_name2<-paste("back_m", ch_available[i], sep="")
              assign(back_m_name2, back_ch_prior[[i]])
              #next # continue the loop to the nect channel

            }

          }

          # post only condition
          if(!exists(paste("back_m_prior", ch_available[i], sep="")) & exists(paste("back_m_post",ch_available[i], sep=""))){

            if(back_ch_post_names[[i]]==paste("back_m_post", ch_available[i], sep="")){
              message("matching background Channels")

              back_m_name2<-paste("back_m", ch_available[i], sep="")
              assign(back_m_name2, back_ch_post[[i]])
              #next # continue the loop to the next channel

            }

          }

          if(exists(paste("back_m_prior", ch_available[i], sep="")) & exists(paste("back_m_post",ch_available[i], sep=""))){

            if(back_ch_post_names[[i]]==paste("back_m_post", ch_available[i], sep="")){
              # message("matching background Channels")
              prior_post_mean<-(back_ch_prior[[i]] + back_ch_post[[i]]) / 2

              back_m_name2<-paste("back_m", ch_available[i], sep="")
              assign(back_m_name2,  prior_post_mean)
              #next # continue the loop to the nect channel

            }

          } # closes 'prior and post available' condition statement

        } # closes the loop -  here have

      }else{ #match_background_Ch=TRUE switch to FALSE

        # 1. prior file only
        # 2. post file only
        # 3. prior and post

        if(is.na(background_post) | !is.na(background_prior)){

          prior_mean<-sum(as.numeric(back_ch_prior)) / (length (unique(back_prior$Ch)))
          back_m<-prior_mean

        }

        if(!is.na(background_post) | is.na(background_prior)){
          post_mean<-sum(as.numeric(back_ch_post)) / (length (unique(back_post$Ch)))
          back_m<-post_mean
        }

        if(!is.na(background_post) & !is.na(background_prior)){

          prior_mean<-sum(as.numeric(back_ch_prior)) / (length (unique(back_prior$Ch)))
          post_mean<-sum(as.numeric(back_ch_post)) / (length (unique(back_post$Ch)))
          back_m<-(prior_mean+post_mean) / 2

        }

      } # end of match background == FALSE (the else part of if statement)
    }# end of background_linear_gr == FALSE

  }# the end of getting the background slopes

  # if match_background_Ch=FALSE then correct all SMR values using back_m (the total avrg)
  # if match_background_Ch=TRUE then correct all SMR values using Chanel specific back_mCh, where Ch is the number of the channel (the total avrg)
  # END -- >>> background

  # if(length(data.SDA)==1){
  #   if(file.exists(data.SDA) | file.exists(paste("./SDA/csv_input_files/", data.SDA[1], sep=""))){ # after running through RMRrepeat - this will be saved in csv input files
  #   	  if(file.exists(paste("./SDA/csv_input_files/", data.SDA[1], sep=""))){
  #         d_SMR<-read.csv(paste("./SDA/csv_input_files/", data.SDA[1], sep=""))
  #       }
  #       if(file.exists(data.SDA)){
  #         d_SMR<-read.csv(data.SDA)
  #       }
  #   	}else{
  #       stop_function<-TRUE
  #       if(stop_function){
  #         stop("Cannot locate the indicated data.SDA data file.")
  #       }
  #   }
  # }

  # SMR, mo2 overnight -----
  # drop any unwanted channels
  if(!is.null(drop_ch[1])){
    n_ch_drop<-length(drop_ch)
    for(i in 1:n_ch_drop){
      d_SMR<-d_SMR[(substr(as.character(d_SMR$Ch), start=3, stop=3)!=drop_ch[i]),]
    }
  }
  d_SMR$Ch<-factor(d_SMR$Ch)

  if(!colnames(d_SMR)[11]=="type"){
    d_SMR$type="SMR"
    d_SMR<-d_SMR[,c("time_frame", "min_start", "r2", "b", "m", "t_min", "t_max", "t_mean" ,"Ch", "DateTime_start", "type", "n_min", "ID_code" )]
  }


  ## choose what to keep and what not for the SMR
  # 2) keep only > 2 min sections for SMR calculations any type  selected above (SMR, pre-, post-shallow slopes) and exclude "SMR-cut", "SMR-cut1", "SMR-cut2"
	d_SMR.1<-d_SMR[d_SMR$n_min>=1,]
	if (nrow(d_SMR.1)==0){
	  warning("All SMR measurements shorter than 1 min !!")
   }else{
    d_SMR<-d_SMR[d_SMR$n_min>=1,]
    d_SMR.2<-d_SMR[d_SMR$n_min<=1,]
      if(nrow(d_SMR.2)==0){
          # message("All SMR/RMR measurements are > 1 min")
        }else{
          message(paste("Using only SMR measurements > 1 min: ", "Not using", nrow(d_SMR.2), " cycles", sep=""))
      }
	}


  # 3) keep only sections with cycles above a set threshold of R2
  if(nrow(d_SMR[d_SMR$r2>=r2_threshold_smr,])<1){
    message(paste("All SMR slope R2 are below the specified threshold. Lower the threshold and rerun the function. NOTE that the lowest R2 is ", min(d_SMR$r2), sep=""))
  }else{
    d_SMR<-d_SMR[d_SMR$r2>=r2_threshold_smr,]
  }

  # first get MO2 values in kg
  d_SMR$bw<-NA
  d_SMR$mo2<-NA
  d_SMR$ID<-NA
  d_SMR$resp.V<-NA
  d_SMR$background_slope<-NA
  d_SMR$scaling_exponent<-scaling_exponent_smr

  for(i in 1:4){
    if(any(grepl(as.character(i),as.character(d_SMR$Ch)))){
      bw.val<-BW.animal[i]
      ID<-AnimalID[i]
      nameCh<-paste("Ch",i,sep="")
      resp.Vol<-resp.V[i]
      n.row<-which(d_SMR$Ch==nameCh)
      d_SMR$bw[n.row]<-bw.val
      d_SMR$ID[n.row]<-as.character(ID)
      d_SMR$resp.V[n.row]<-resp.Vol

    }
  }


  # START -- >>> background corrections SMR

  # MO2 values MR in mgO2/min/kg - background corrected
  # ------- ###
  # !!! ONLY FOR SMR background corrections
  # 1. if background files (either prior or post, or both) are provided we account for it
  # 1.1 if match_background_Ch=FALSE then use back_m to correct each mo2 value
  # 1.2 if match_background_Ch=TRUE then use back_m[Ch] to correct each mo2 value for each channel
  #  2. if background files are NOT provided we DONT acount for any background respiration

  if(!is.na(strptime(d_SMR$DateTime_start[1], format = date_format[1], tz = ""))){
    d_SMR$DateTime_start<- strptime(d_SMR$DateTime_start, format = date_format, tz = "")
  }
  if(!is.na(strptime(d_SMR$DateTime_start[1], format = date_format[2], tz = ""))){
    d_SMR$DateTime_start<- strptime(d_SMR$DateTime_start, format = date_format[2], tz = "")
  }

  # 000 Jan 4 addition - account for background using linear regression
  # making predictions
  if(!is.null(background_gr) & match_background_Ch==FALSE){
    background_slopes<-data.frame(DateTime_start = as.POSIXct(d_SMR$DateTime_start))

    if(background_gr == "linear"){
      background_slopes$back_m<-predict(back_regression, background_slopes)
    }
    if(background_gr == "exp"){
      background_slopes$back_m<-exp(predict(back_regression, background_slopes))
    }

    if (local_path | !dir.exists("./SDA/plots_background")){
      plotname.backgr<-paste( filename.SMR,"_PLOT_BACKGROUND.png", sep="")
    }else{
      plotname.backgr<-paste("./SDA/plots_background/", filename.SMR,"_PLOT_BACKGROUND.png", sep="")
    }

    for (i in 1:nrow(d_SMR)){
      d_SMR$mo2[i]<-(( d_SMR$m[i]*(d_SMR$resp.V[i]-d_SMR$bw[i])) - (background_slopes$back_m[i] * d_SMR$resp.V[i])) /(d_SMR$bw[i])#^scaling_exponent_smr) # units mgO2 kg-1 min-1 -CORRECTED for back resp
    }

    background_slopes[,3:4]<-d_SMR[,c(5,9)]

    background_slope_plot<-ggplot(data=background_slopes, aes(x=DateTime_start, m, colour=Ch))+
      geom_point(size=2, pch=21)+
      geom_point(aes(x=DateTime_start, back_m), colour="black", fill="grey", alpha=0.5, pch=19, size=1)+
      theme_bw()+
      facet_grid(Ch~.)

    # plotname.backgr<-paste( filename.SMR,"_PLOT_BACKGROUND.png", sep="")
    png(plotname.backgr, width=4, height=8, res=300, units="in")
      print(background_slope_plot)
    dev.off()

    d_SMR$background_slope <- background_slopes$back_m

  }

  if(!is.null(background_gr) & match_background_Ch==TRUE){

    if (local_path | !dir.exists("./SDA/plots_background")){
      plotname.backgr<-paste( filename.SMR,"_PLOT_BACKGROUND.png", sep="")
    }else{
      plotname.backgr<-paste("./SDA/plots_background/", filename.SMR,"_PLOT_BACKGROUND.png", sep="")
    }

    message("SMR corrected for background: using Ch specific average background")
    background_slopes<-data.frame(DateTime_start = d_SMR$DateTime_start)

    for (i in 1:nrow(d_SMR)){

      # print(c(d_SMR$bw[1],d_SMR$Ch[i]))
      if(substr(d_SMR$Ch[i], start=3, stop=3) == "1"){
  	      if(background_gr == "linear"){
    	      background_slopes$back_m[i]<-predict(back_regression1, data.frame(DateTime_start = d_SMR$DateTime_start[i]))
            back_m1<-predict(back_regression1, data.frame(DateTime_start = d_SMR$DateTime_start[i]))
  	      }
  	      if(background_gr == "exp"){
    	      background_slopes$back_m[i]<-exp(predict(back_regression1, data.frame(DateTime_start = d_SMR$DateTime_start[i])))
            back_m1<-exp(predict(back_regression1, data.frame(DateTime_start = d_SMR$DateTime_start[i])))
  	      }
        d_SMR$mo2[i]<-(( d_SMR$m[i]*(d_SMR$resp.V[i]-d_SMR$bw[i])) - (back_m1 * d_SMR$resp.V[i])) /(d_SMR$bw[i])#^scaling_exponent_smr) # units mgO2 kg-1 min-1 -CORRECTED for back resp
        # message("Channel: 1,  bacterial respiration slope: ", back_m1)
      }

      if(substr(d_SMR$Ch[i], start=3, stop=3) == "2"){
  	     if(background_gr == "linear"){
    	      background_slopes$back_m[i]<-predict(back_regression2, data.frame(DateTime_start = d_SMR$DateTime_start[i]))
            back_m2<-predict(back_regression2, data.frame(DateTime_start = d_SMR$DateTime_start[i]))
  	      }
  	      if(background_gr == "exp"){
    	      background_slopes$back_m[i]<-exp(predict(back_regression2, data.frame(DateTime_start = d_SMR$DateTime_start[i])))
            back_m2<-exp(predict(back_regression2, data.frame(DateTime_start = d_SMR$DateTime_start[i])))
  	      }
        d_SMR$mo2[i]<-(( d_SMR$m[i]*(d_SMR$resp.V[i]-d_SMR$bw[i])) - (back_m2 * d_SMR$resp.V[i])) /(d_SMR$bw[i])#^scaling_exponent_smr) # units mgO2 kg-1 min-1 -CORRECTED for back resp
        # message("Channel: 2,  bacterial respiration slope: ", back_m2)
      }

      if(substr(d_SMR$Ch[i], start=3, stop=3) == "3"){
  	     if(background_gr == "linear"){
    	      background_slopes$back_m[i]<-predict(back_regression3, data.frame(DateTime_start = d_SMR$DateTime_start[i]))
            back_m3<-predict(back_regression3, data.frame(DateTime_start = d_SMR$DateTime_start[i]))
  	      }
  	      if(background_gr == "exp"){
    	      background_slopes$back_m[i]<-exp(predict(back_regression3, data.frame(DateTime_start = d_SMR$DateTime_start[i])))
            back_m3<-exp(predict(back_regression3, data.frame(DateTime_start = d_SMR$DateTime_start[i])))
  	      }
        d_SMR$mo2[i]<-(( d_SMR$m[i]*(d_SMR$resp.V[i]-d_SMR$bw[i])) - (back_m3 * d_SMR$resp.V[i])) /(d_SMR$bw[i])#^scaling_exponent_smr) # units mgO2 kg-1 min-1 -CORRECTED for back resp
        # message("Channel: 3,  bacterial respiration slope: ", back_m3)
      }

      if(substr(d_SMR$Ch[i], start=3, stop=3) == "4"){
  	     if(background_gr == "linear"){
    	      background_slopes$back_m[i]<-predict(back_regression4, data.frame(DateTime_start = d_SMR$DateTime_start[i]))
            back_m4<-predict(back_regression4, data.frame(DateTime_start = d_SMR$DateTime_start[i]))
  	      }
  	      if(background_gr == "exp"){
    	      background_slopes$back_m[i]<-exp(predict(back_regression4, data.frame(DateTime_start = d_SMR$DateTime_start[i])))
            back_m4<-exp(predict(back_regression4, data.frame(DateTime_start = d_SMR$DateTime_start[i])))
  	      }
        d_SMR$mo2[i]<-(( d_SMR$m[i]*(d_SMR$resp.V[i]-d_SMR$bw[i])) - (back_m4 * d_SMR$resp.V[i])) /(d_SMR$bw[i])#^scaling_exponent_smr) # units mgO2 kg-1 min-1 -CORRECTED for back resp
        # message("Channel: 4,  bacterial respiration slope: ", back_m4)
      }
    }# end of the loop

    background_slopes[,3:4]<-d_SMR[,c(5,9)]
    background_slope_plot<-ggplot(data=background_slopes, aes(x=DateTime_start, m, colour=Ch))+
      geom_point(size=2, pch=21)+
      geom_point(aes(x=DateTime_start, back_m), colour="black", fill="grey", alpha=0.5, pch=19, size=1)+
      theme_bw()+
      facet_grid(Ch~.)

    # WORK - SAVE this plot
    # plotname.backgr<-paste( filename.SMR,"_PLOT_BACKGROUND.png", sep="")
    png(plotname.backgr, width=4, height=8, res=300, units="in")
      print(background_slope_plot)
    dev.off()

    d_SMR$background_slope<- background_slopes$back_m
  }

  # 1.1 if background files (either prior or post, or both) are provided and its one overall mean value (back_m)
  if ((( !is.null(background_post) | !is.null(background_prior)) & match_background_Ch==FALSE) & is.null(background_slope) & is.null(background_gr)){
    message("SMR corrected for background: used a mean (prior and/or post) background measurements | mean bacterial respiration slope: ", back_m, " ~" , round((back_m*100)/(mean(d_SMR[d_SMR$m <= quantile(d_SMR$m, 0.5, na.rm=TRUE), "m"], na.rm=TRUE)), 2), " %")

    for (i in 1:nrow(d_SMR)){
      d_SMR$mo2[i]<-(( d_SMR$m[i]*(d_SMR$resp.V[i]-d_SMR$bw[i])) - (back_m * d_SMR$resp.V[i])) /(d_SMR$bw[i])#^scaling_exponent_smr) # units mgO2 kg-1 min-1 -CORRECTED for back resp
    }

    d_SMR$background_slope<- back_m
  }

  ## if background slope and volume are specifically provided, then use those! this alos overrides the background prior and post argument.
  # all channels with the same slope
  if (!is.null(background_slope)){
    message("SMR corrected for background: used a common manually provided background slope for all channels | bacterial respiration slope: ", background_slope, " ~" , round((background_slope*100)/(mean(d_SMR[d_SMR$m <= quantile(d_SMR$m, 0.5, na.rm=TRUE), "m"], na.rm=TRUE)), 2), " %")
    for (i in 1:nrow(d_SMR)){
      d_SMR$mo2[i]<-(( d_SMR$m[i]*(d_SMR$resp.V[i]-d_SMR$bw[i])) - (background_slope * background.V)) /(d_SMR$bw[i])#^scaling_exponent_smr) # units mgO2 kg-1 min-1 -CORRECTED for back resp
    }
    d_SMR$background_slope<- paste(round(background_slope,2), "_BackVol=", background.V, sep="")
  }

  # 1.2 if background files are provided and its channel specific
  if ((!is.null(background_post) | !is.null(background_prior)) & match_background_Ch==TRUE & is.null(background_gr)){
    message("SMR corrected for background: using Ch specific average background")

    for (i in 1:nrow(d_SMR)){

      # print(c(d_SMR$bw[1],d_SMR$Ch[i]))
      if(substr(d_SMR$Ch[i], start=3, stop=3) == "1"){
        d_SMR$mo2[i]<-(( d_SMR$m[i]*(d_SMR$resp.V[i]-d_SMR$bw[i])) - (back_m1 * d_SMR$resp.V[i])) /(d_SMR$bw[i])#^scaling_exponent_smr) # units mgO2 kg-1 min-1 -CORRECTED for back resp
        # message("Channel: 1,  bacterial respiration slope: ", back_m1)
        d_SMR$background_slope[i]<-back_m1
      }

      if(substr(d_SMR$Ch[i], start=3, stop=3) == "2"){
        d_SMR$mo2[i]<-(( d_SMR$m[i]*(d_SMR$resp.V[i]-d_SMR$bw[i])) - (back_m2 * d_SMR$resp.V[i])) /(d_SMR$bw[i])#^scaling_exponent_smr) # units mgO2 kg-1 min-1 -CORRECTED for back resp
        # message("Channel: 2,  bacterial respiration slope: ", back_m2)
        d_SMR$background_slope[i]<-back_m2
      }

      if(substr(d_SMR$Ch[i], start=3, stop=3) == "3"){
        d_SMR$mo2[i]<-(( d_SMR$m[i]*(d_SMR$resp.V[i]-d_SMR$bw[i])) - (back_m3 * d_SMR$resp.V[i])) /(d_SMR$bw[i])#^scaling_exponent_smr) # units mgO2 kg-1 min-1 -CORRECTED for back resp
        # message("Channel: 3,  bacterial respiration slope: ", back_m3)
        d_SMR$background_slope[i]<-back_m3
      }

      if(substr(d_SMR$Ch[i], start=3, stop=3) == "4"){
        d_SMR$mo2[i]<-(( d_SMR$m[i]*(d_SMR$resp.V[i]-d_SMR$bw[i])) - (back_m4 * d_SMR$resp.V[i])) /(d_SMR$bw[i])#^scaling_exponent_smr) # units mgO2 kg-1 min-1 -CORRECTED for back resp
        # message("Channel: 4,  bacterial respiration slope: ", back_m4)
        d_SMR$background_slope[i]<-back_m4
      }
    }# end of the loop

  }

  # 2. if background files are not provided
  if ((is.null(background_post) & is.null(background_prior)) & is.null(background_slope)){
    message("SMR/RMR: NO correction for background respiration")

    for (i in 1:nrow(d_SMR)){
      d_SMR$mo2[i]<-d_SMR$m[i]*(d_SMR$resp.V[i]-d_SMR$bw[i])/(d_SMR$bw[i])#^scaling_exponent_smr) # units mgO2 kg-1 min-1
    }
  }
  ## end of SMR corrections for body size and background

  # END -- >>> background corrections SMR

  # find SMR estimates using all methods provided by Chabot
  # --- #
  # 1. frequency plot

  p_freq<-ggplot(d_SMR, aes(x=mo2))+
    geom_histogram(bins = 60 , color = "black", fill = "gray") +
    facet_grid(Ch~.)+
    theme_classic()+
    ggtitle(data.SDA[1])

  d_SMR$DateTime_start<-as.character(d_SMR$DateTime_start)
  # the lowest 10 values after removal of 5 lowest
  a0<-d_SMR[,2:ncol(d_SMR)] %>%
    group_by(Ch)%>%
    top_n(-15, mo2)%>%
    arrange(Ch,mo2)

  # the 5 lowest/ exluded
  a00<-a0 %>%
    group_by(Ch)%>%
    top_n(-5, mo2)%>%
    arrange(Ch,mo2)

  # the 10 lowest/ exluding the 5 lowest in the df
  a<-a0 %>%
    group_by(Ch)%>%
    top_n(10, mo2)%>%
    arrange(Ch,mo2)

  min10_MO2<-as.data.frame(a)

  min10_mean<-min10_MO2%>%
    group_by(Ch)%>%
    summarize(mean_temp=mean(t_mean), sd_temp=sd(t_mean), mean_mo2=mean(mo2), sd_mo2=sd(mo2),
              cv_mo2 = sd(mo2)/(sqrt(10)), n=length(mo2))

  min10_plot<-ggplot(data=d_SMR, aes(x=min_start, y=mo2))+
    geom_point(size=1)+
    geom_point(data=a00, aes(x=min_start, y=mo2), color="red", pch=19, size=3)+
    geom_point(data=min10_MO2, aes(x=min_start, y=mo2), colour="green4",size=3, alpha=0.7)+
    geom_line(linewidth=0.5, alpha=0.)+
		theme_light()+
    ggtitle("lowest 10 VALUES  excluding 5 smallest")+
    theme(legend.position="top")+
    facet_grid(Ch~.)

  # the 10% percentile
  # perc=c(0.1,0.15, 0.2)

  a2<-as.data.frame(matrix(ncol=3, nrow=0))

  if(any(is.na(d_SMR$mo2) | is.nan(d_SMR$mo2))){
    stop_function<-TRUE
    if(stop_function) {
      print(min10_plot)
      stop("See exported plot: some MO2 values are NA or NaN, cannot estimate SMR")
    }
  }

  for (i in unique(d_SMR$Ch)){
    split_temp<-as.data.frame(split(d_SMR, d_SMR$Ch)[i])
    colnames(split_temp)<-colnames(d_SMR)
    quant10<-quantile(split_temp$mo2, 0.1, na.rm=T)
    quant15<-quantile(split_temp$mo2, 0.15, na.rm=T)
    quant20<-quantile(split_temp$mo2, 0.2, na.rm=T)

    a2<-rbind(a2,as.data.frame(t(c(as.character(split_temp$Ch[1]), "10%", as.numeric(quant10)))))
    a2<-rbind(a2,as.data.frame(t(c(as.character(split_temp$Ch[1]), "15%", as.numeric(quant15)))))
    a2<-rbind(a2,as.data.frame(t(c(as.character(split_temp$Ch[1]), "20%", as.numeric(quant20)))))
  }

  colnames(a2)<-c("Ch","quantiles", "mo2_perc")
  a2$mo2_perc<-as.numeric(as.character( a2$mo2_perc))

  quantile_smr<-spread(a2, key=quantiles, value=mo2_perc)
  a2$quantiles<-as.factor(a2$quantiles)

  row_ch1_10perc<-which(a2$Ch=="Ch1" & a2$quantiles=="10%")
  row_ch2_10perc<-which(a2$Ch=="Ch2" & a2$quantiles=="10%")
  row_ch3_10perc<-which(a2$Ch=="Ch3" & a2$quantiles=="10%")
  row_ch4_10perc<-which(a2$Ch=="Ch4" & a2$quantiles=="10%")

  row_ch1_15perc<-which(a2$Ch=="Ch1" & a2$quantiles=="15%")
  row_ch2_15perc<-which(a2$Ch=="Ch2" & a2$quantiles=="15%")
  row_ch3_15perc<-which(a2$Ch=="Ch3" & a2$quantiles=="15%")
  row_ch4_15perc<-which(a2$Ch=="Ch4" & a2$quantiles=="15%")

  row_ch1_20perc<-which(a2$Ch=="Ch1" & a2$quantiles=="20%")
  row_ch2_20perc<-which(a2$Ch=="Ch2" & a2$quantiles=="20%")
  row_ch3_20perc<-which(a2$Ch=="Ch3" & a2$quantiles=="20%")
  row_ch4_20perc<-which(a2$Ch=="Ch4" & a2$quantiles=="20%")

  mo2_lab<-bquote(MO[2]~(mgO[2]~min^-1~ kg^-1))

  min_percPlot<-ggplot(data=d_SMR, aes(x=min_start, y=mo2))+
			geom_point(size=1)+
		  geom_point(data=d_SMR[which(d_SMR$Ch=="Ch1" & d_SMR$mo2<=a2$mo2_perc[row_ch1_20perc]),], aes(x=min_start, y=mo2), colour="purple",size=3)+
		  geom_point(data=d_SMR[which(d_SMR$Ch=="Ch2" & d_SMR$mo2<=a2$mo2_perc[row_ch2_20perc]),], aes(x=min_start, y=mo2), colour="purple",size=3)+
		  geom_point(data=d_SMR[which(d_SMR$Ch=="Ch3" & d_SMR$mo2<=a2$mo2_perc[row_ch3_20perc]),], aes(x=min_start, y=mo2), colour="purple",size=3)+
		  geom_point(data=d_SMR[which(d_SMR$Ch=="Ch4" & d_SMR$mo2<=a2$mo2_perc[row_ch4_20perc]),], aes(x=min_start, y=mo2), colour="purple",size=3)+
		  geom_point(data=d_SMR[which(d_SMR$Ch=="Ch1" & d_SMR$mo2<=a2$mo2_perc[row_ch1_15perc]),], aes(x=min_start, y=mo2), colour="blue",size=3)+
		  geom_point(data=d_SMR[which(d_SMR$Ch=="Ch2" & d_SMR$mo2<=a2$mo2_perc[row_ch2_15perc]),], aes(x=min_start, y=mo2), colour="blue",size=3)+
		  geom_point(data=d_SMR[which(d_SMR$Ch=="Ch3" & d_SMR$mo2<=a2$mo2_perc[row_ch3_15perc]),], aes(x=min_start, y=mo2), colour="blue",size=3)+
		  geom_point(data=d_SMR[which(d_SMR$Ch=="Ch4" & d_SMR$mo2<=a2$mo2_perc[row_ch4_15perc]),], aes(x=min_start, y=mo2), colour="blue",size=3)+
		  geom_point(data=d_SMR[which(d_SMR$Ch=="Ch1" & d_SMR$mo2<=a2$mo2_perc[row_ch1_10perc]),], aes(x=min_start, y=mo2), colour="darkturquoise",size=3)+
		  geom_point(data=d_SMR[which(d_SMR$Ch=="Ch2" & d_SMR$mo2<=a2$mo2_perc[row_ch2_10perc]),], aes(x=min_start, y=mo2), colour="darkturquoise",size=3)+
		  geom_point(data=d_SMR[which(d_SMR$Ch=="Ch3" & d_SMR$mo2<=a2$mo2_perc[row_ch3_10perc]),], aes(x=min_start, y=mo2), colour="darkturquoise",size=3)+
		  geom_point(data=d_SMR[which(d_SMR$Ch=="Ch4" & d_SMR$mo2<=a2$mo2_perc[row_ch4_10perc]),], aes(x=min_start, y=mo2), colour="darkturquoise",size=3)+
		  # geom_p# geom_line(size=0.5, alpha=0.5)+
			theme_light()+
			# ggtitle("10th percentile")+
			theme(legend.position="top")+
		  ylab(mo2_lab)+
		  xlab("Time in trial (min)")+
			facet_grid(Ch~.)

		png(plotname.freq, width = 10, height = 10, units="in", res=200)
			grid.arrange( min10_plot, min_percPlot, ncol=1, nrow=2)
		dev.off()

  # MLND algorithm from suppl JFB issue 88 Chabot, Steffensen, and Farrell 2016
  # SMR data frames split based on the channel
  if (length(unique(d_SMR$Ch))>1){
    Ch.data.smr<-split(d_SMR, d_SMR$Ch)
  }else{
    Ch.data.smr<-1
  }

  for (i in 1:length(unique(d_SMR$Ch))){

    if (length(unique(d_SMR$Ch))==1){
      Y0<-d_SMR
    }else{
      Y0<-as.data.frame(Ch.data.smr[i])
    }

    colnames(Y0)<-c("time_frame","min_start", "r2" ,"b", "m" , "t_min", "t_max", "t_mean", "Ch", "DateTime_start","mo2_type","n_min", "ID_code","bw", "mo2","ID")

    mo2<-Y0$mo2
    BW<-Y0$bw[1]
    ID<-Y0$ID[1]
    t_min<-min(Y0$t_min) # min of min temp values
    t_max<-max(Y0$t_max) # max of max temp values
    t_mean<-mean(Y0$t_mean) # mean of mean temp values
    N_mo2<-length(mo2)

    Y.Ch<-as.character(Y0$Ch[1])
		if (local_path | !dir.exists("SDA")){
      plotname.mlnd<-paste( filename.SMR,"_", Y.Ch, "_MLND_SMR_analysed.png", sep="")
  	}else{
      plotname.mlnd<-paste("./SDA/plots_mlnd_SMR/", filename.SMR,"_", Y.Ch, "_MLND_SMR_analysed.png", sep="")
  	}


    # START -->>> MLND
    # loop for MLND function
    if(MLND){

      mlnd.m1<-Mclust(mo2, G=1:4, verbose = verbose.MLND)
      m1<-densityMclust(mo2, verbose = verbose.MLND) # same as mlnd
      clas.m1<-m1$classification #
      clas.m1.t <- as.data.frame(table(clas.m1))
      clas.m1.t$clas.m1 <- as.numeric(levels(clas.m1.t$clas.m1))
      boot1 <- MclustBootstrap(m1, nboot = 1000, type = "bs", verbose = verbose.MLND)


      valid <- clas.m1.t$Freq>=0.1*length(mo2)   # assuming that the lower class contains 10% of the MO2 values recorded

      if(max(clas.m1)==1 ){

        message("ONLY ONE MLND")
        if (valid[1]){
          message("MLND class 1 = valid")
          valid.clas<-1
        }
        png(plotname.mlnd, width = 12, height = 9, units="in", res=200)
        par(mfrow=c(3,4))
        #~ 					plot(m1, what= "BIC")
        plot(mlnd.m1, what="classification")
        plot(m1, what = "density", data = mo2, breaks = 30)
        plot(m1, what = "diagnostic", type = "cdf")
        plot(m1, what = "diagnostic", type = "qq")
        plot(Y0$mo2~Y0$min_start, cex=0.5)
        points(Y0$mo2[which(clas.m1==1)]~Y0$min_start[which(clas.m1==1)], col="red", cex=1)
        #~ 				plot(boot1, what = "mean", show.confint=TRUE)
        plot(1, type="n", axes=F, xlab="", ylab="");
        text(1, 1, Y.Ch ,cex = 3)
        #~ 				plot(boot1, what = "pro")
        dev.off()

      }else{

        if (!is.na(boot1$variance[1])){
          sum.boot1<-summary(boot1, what="ci")
        }

        if(max(clas.m1)==2){
          if (valid[1]){
            message("MLND class 1 = valid")
            valid.clas<-1
            CIup.mlnd<-sum.boot1$mean[2]
            CIlo.mlnd<-sum.boot1$mean[1]

          }else{
            valid.clas<- min(clas.m1.t$clas.m1[valid])
            message(paste("MLND lowest valid class = ", valid.clas, sep=""))

            if(valid.clas==2){
              CIup.mlnd<-sum.boot1$mean[4]
              CIlo.mlnd<-sum.boot1$mean[3]
            }
          }

          png(plotname.mlnd, width = 12, height = 9, units="in", res=200)
          par(mfrow=c(3,4))
          #~ 					plot(m1, what= "BIC")
          plot(mlnd.m1, what="classification")
          plot(m1, what = "density", data = mo2, breaks = 30)
          plot(m1, what = "diagnostic", type = "cdf")
          plot(m1, what = "diagnostic", type = "qq")
          plot(Y0$mo2~Y0$min_start, cex=0.5)
          points(Y0$mo2[which(clas.m1==1)]~Y0$min_start[which(clas.m1==1)], col="red", cex=1)
          plot(boot1, what = "mean", show.confint=TRUE)
          plot(1, type="n", axes=F, xlab="", ylab="");
          text(1, 1, "" ,cex = 3)
          plot(1, type="n", axes=F, xlab="", ylab="");
          text(1, 1, Y.Ch ,cex = 3)
          plot(boot1, what = "pro")
          dev.off()

        }

        if(max(clas.m1)>=3){
          if (valid[1]){
            message("MLND class 1 = valid")
            valid.clas<-1
            CIup.mlnd<-sum.boot1$mean[2]
            CIlo.mlnd<-sum.boot1$mean[1]

          }else{
            valid.clas<- min(clas.m1.t$clas.m1[valid])
            message(paste("MLND lowest valid class = ", valid.clas, sep=""))

            if(valid.clas==2){
              CIup.mlnd<-sum.boot1$mean[4]
              CIlo.mlnd<-sum.boot1$mean[3]
            }else{
              CIup.mlnd<-sum.boot1$mean[6]
              CIlo.mlnd<-sum.boot1$mean[5]
            }
          }


          png(plotname.mlnd, width = 12, height = 12, units="in", res=200)
          par(mfrow=c(4,4))
          #~ 					plot(m1, what= "BIC")
          plot(mlnd.m1, what="classification")
          plot(m1, what = "density", data = mo2, breaks = 30)

          if (!max(clas.m1)>3){
            plot(m1, what = "diagnostic", type = "cdf")
            plot(m1, what = "diagnostic", type = "qq")
          }else{
            message("More than 3 MLND classes, no diagnostic plots")
          }

          plot(Y0$mo2~Y0$min_start, cex=0.5)
          points(Y0$mo2[which(clas.m1==1)]~Y0$min_start[which(clas.m1==1)], col="red", cex=1)
          plot(boot1, what = "mean", show.confint=TRUE)
          plot(1, type="n", axes=F, xlab="", ylab="");
          text(1, 1, Y.Ch ,cex = 3)
          plot(boot1, what = "pro")
          dev.off()

        }

      }

      distr <- mo2[m1$classification==valid.clas]
      mlnd <- m1$parameters$mean[valid.clas]
      CVmlnd <- sd(distr)/mlnd * 100
      Nmlnd<- length(distr)



    }else{ # end of if(MLND=TRUE){
      if(substr(Y.Ch, start=3, stop=3)=="1"){
        message("MLND not calculated")
      }
      distr <- NA
      mlnd <- 0
      CVmlnd <- 0
      Nmlnd <- 0

    }
    # END -- >>> MNLD
    # https://cran.r-project.org/web/packages/mclust/vignettes/mclust.html#classification
    # mclust: R package for model-based clustering, classification, and density estimation based on finite normal mixture modelling. It provides functions for parameter estimation via the EM algorithm for normal mixture models with a variety of covariance structures, and functions for simulation from these models.
    # read about EM http://www.statisticshowto.com/em-algorithm-expectation-maximization/




    values_smr<-as.data.frame(t(c(filename.SMR, ID, Y.Ch, BW, t_min, t_max, t_mean, N_mo2,
                                  min10_mean$mean_mo2[i], min10_mean$sd_mo2[i], min10_mean$cv_mo2[i],
                                  quantile_smr[i,2],  quantile_smr[i,3], quantile_smr[i,4],
                                  mlnd, CVmlnd, Nmlnd)))

    colnames(values_smr)<-c("filename", "ID", "Ch", "BW","t_min","t_max", "t_mean", "N_mo2", #8
                            "smr_mean10minVal","smr_SD10minVal", "smr_CV10minVal", "SMR_low10quant","SMR_low15quant","SMR_low20quant", #6
                            "smr_mlnd", "smr_CVmlnd", "smr_Nmlnd")#1

    newdata.smr<-rbind(newdata.smr, values_smr)

  }# end of MNLND calculatons. Loops throgh every channel

  pd<-as.data.frame(newdata.smr[,c(3,9, 12,13,14,15, 8, 17)])

  df.1 <- data.frame(Ch=unlist(pd[,1], use.names=FALSE), smr_val=unlist(pd[,2], use.names=FALSE), smr_method="smr_mean10minVal", N_mo2=unlist(pd[,7], use.names=FALSE))
  df.2 <- data.frame(Ch=unlist(pd[,1], use.names=FALSE), smr_val=unlist(pd[,3], use.names=FALSE), smr_method="SMR_low10quant", N_mo2=unlist(pd[,7], use.names=FALSE))
  df.3 <- data.frame(Ch=unlist(pd[,1], use.names=FALSE), smr_val=unlist(pd[,4], use.names=FALSE), smr_method="SMR_low15quant", N_mo2=unlist(pd[,7], use.names=FALSE))
  df.4 <- data.frame(Ch=unlist(pd[,1], use.names=FALSE), smr_val=unlist(pd[,5], use.names=FALSE), smr_method="SMR_low20quant", N_mo2=unlist(pd[,7], use.names=FALSE))
  df.5 <- data.frame(Ch=unlist(pd[,1], use.names=FALSE), smr_val=unlist(pd[,6], use.names=FALSE), smr_method="smr_mlnd", N_mo2=unlist(pd[,8], use.names=FALSE))
  plot_d<-rbind(df.1, df.2, df.3, df.4, df.5 )

  plot_d$smr_val<-as.numeric(as.character(plot_d$smr_val))

  smr_meth_p<-
    ggplot(plot_d, aes(y=smr_val, x=factor(smr_method), colour=factor(Ch), group=factor(Ch), label=as.character(N_mo2)))+
    geom_point(size=7,pch=21)+
    geom_text(color="black")+
    geom_line()+
    theme_classic()+
    ylab(expression(SMR~(mg~O[2]~kg^-1~min^-1)))+
    xlab("SMR test ")+
    theme(axis.text.y=element_text(size=20, colour= 'black'),
          axis.text.x=element_text(size=15, colour= 'black', angle=90, hjust=1),
          axis.line.y=element_line(colour = 'black',size=0.5),
          axis.line.x=element_line(colour = 'black',size=0.5),
          axis.ticks.y=element_line(size=0.5),
          axis.ticks.x=element_line(size=0))+
    theme(axis.title.y=element_text(size=15),
          axis.title.x=element_text(size=15),
          panel.border = element_rect(size=0.9,linetype = "solid",fill=NA, colour = "black"))

  png(plotname.smr.meth, width=6, height=5, units="in", res=200)
    print(smr_meth_p)
  dev.off()


  # save data:
  if (local_path | !dir.exists("SDA")){
    filename.smr<-paste( gsub('.{12}$', '', data.SDA[1]), "SMR_analyzed.csv", sep='')
    filename.MR<-paste( gsub('.{12}$', '', data.SDA[1]), "MR_analyzed.csv", sep='')
  }else{
    filename.smr<-paste("./SDA/csv_analyzed_SMR/",gsub('.{12}$', '', data.SDA[1]), "SMR_analyzed.csv", sep='')
    filename.MR<-paste("./SDA/csv_analyzed_MR/", gsub('.{12}$', '', data.SDA[1]), "MR_analyzed.csv", sep='')
  }

  lst <- lapply(newdata.smr, unlist)
  newdata.smr <- (data.frame(lapply(lst, `length<-`, max(lengths(lst)))))

  write.csv(file=filename.smr, d_SMR, row.names=FALSE)
  write.csv(file=filename.MR, newdata.smr, row.names=FALSE)

  # message("Save SMR, and MR files")
  # print(head(newdata.smr))
  # print(head(d_SMR))
  # return(newdata.smr)

  # --- >>> end of SMR_clac == TRUE


  if (SMR_calc==FALSE){

    if (SMR_calc==FALSE & is.null(SMR_vals) & is.null(analyzed_MR)){
      stop("SMR analysis: must provide 'analyzed_MR' file or specific 'SMR_vals'")
    }

    if(!is.null(analyzed_MR)){
      if(file.exists(analyzed_MR) | file.exists(paste("./SDA/csv_input_files/", analyzed_MR, sep=""))){ # after running through RMRrepeat - this will be saved in csv input files
      	  if(file.exists(paste("./SDA/csv_input_files/", analyzed_MR, sep=""))){
            newdata.smr<-read.csv(paste("./SDA/csv_input_files/", analyzed_MR, sep=""))
          }
          if(file.exists(analyzed_MR)){
            newdata.smr<-read.csv(analyzed_MR)
          }
      	}else{
          stop_function<-TRUE
          if(stop_function){
            stop("SMR analysis: cannot locate the indicated analyzed_MR data file with SMR values.")
          }
      }
    }

    # newdata.smr<-read.csv(analyzed_MR)
    message("SMR analysis: using estimated SMR values from previosly analyzed provided data")
  }


  if (local_path | !dir.exists("SDA")){
    plotname.sda.data<-	paste(gsub('.{4}$', '', data.SDA[1]),"_SDA_hourly_PLOT.png", sep='')
    SDAdata_name<-paste(gsub('.{4}$', '', data.SDA[1]),"_SDA_analyzed.csv", sep='')
    SDAhrlydata_name<-paste(gsub('.{4}$', '', data.SDA[1]),"_SDA_hrly_analyzed.csv", sep='')
    SDAhrlydata_name_wDELAY<-paste(gsub('.{4}$', '', data.SDA[1]),"_SDA_hrly_wDELAY_analyzed.csv", sep='')

  }else{
    plotname.sda.data<-	paste("./SDA/plots_SDA_hourly/",gsub('.{4}$', '', data.SDA[1]),"_SDA_hourly_PLOT.png", sep='')
    SDAdata_name<-paste("./SDA/csv_analyzed_SDA/", gsub('.{4}$', '', data.SDA[1]),"_SDA_analyzed.csv", sep='')
    SDAhrlydata_name<-paste("./SDA/csv_analyzed_SDA_hrly/", gsub('.{4}$', '', data.SDA[1]),"_SDA_hrly_analyzed.csv", sep='')
    SDAhrlydata_name_wDELAY<-paste("./SDA/csv_analyzed_SDA_hrly/", gsub('.{4}$', '', data.SDA[1]),"_SDA_hrly_wDELAY_analyzed.csv", sep='')
  }

  # SDA
  # 1. get the mean of each hour
  d_SMR$hour<-as.factor(floor(d_SMR$min_start/60))
  d_SMR$ID<-as.factor(as.character(d_SMR$ID))
  d_SMR_original<-d_SMR
  # handling delay for SDA analysis
  # split data frames (full data set)
  if(length(unique(d_SMR$Ch))>1){
    Ch.data.sda.full<-split(d_SMR, d_SMR$Ch)
  }else{
    Ch.data.sda.full<-d_SMR
  }

  for(h in 1:length(Ch.data.sda.full)){
    # as.data.frame(Ch.data.sda.full[[h]])
    h.Ch.data<-data.frame(Ch.data.sda.full[h])
    names(h.Ch.data)<-names(d_SMR)
    h.Ch<-as.character(h.Ch.data$Ch[1])
    h.Ch.data<-h.Ch.data[c(h.Ch.data$min_start/60) >= handling_delay[as.numeric(substr(h.Ch, start=3, stop=3))] , ]

    if(h == 1){
      Ch.data.sda.full_temp<-h.Ch.data
    }else{
      Ch.data.sda.full_temp<-rbind(Ch.data.sda.full_temp, h.Ch.data)
    } # take the experimental data that is passed the handling delay if any

  }

  d_SMR<-Ch.data.sda.full_temp

  d_SMRsum<-d_SMR %>%
    group_by(Ch, ID, resp.V, bw, hour) %>%
    summarize(mo2_mean = mean(mo2, na.rm=TRUE),
              mo2_sd = sd(mo2, na.rm=TRUE),
              mo2_min = min(mo2, na.rm=TRUE),
              t_mean=mean(t_mean),
              n=length (mo2),
              mo2_max = max(mo2, na.rm=TRUE),
              .groups = "keep")

  d_SMRsum$hour<-as.numeric(as.character(d_SMRsum$hour)) + 0.5 # (use avg 0.5 for the mean hour time)
  d_SMRsum$ID<-as.factor(as.character(d_SMRsum$ID))
  d_SMRsum$resp.V<-as.numeric(as.character(d_SMRsum$resp.V))
  d_SMRsum$bw<-as.numeric(as.character(d_SMRsum$bw))

  n_id<-length(unique(levels(d_SMRsum$ID)))

  # get split dataframes (d_SMRsum)
  if (length(unique(d_SMRsum$Ch))>1){
    Ch.data.sda<-split(d_SMRsum, d_SMRsum$Ch)
  }else{
    Ch.data.sda<-d_SMRsum
  }

  # SDA analysis function on hourly data frame
  SDAdata<-matrix(ncol=15, nrow=0)
  colnames(SDAdata)<-c("ID","SMR","spar", "SDA_integrated", "end_SDA_estim_hr", "SMR_intergrated", "peak_SDA", "time_peak_SDA", "percentSMR_peak_SDA", "MO2_SDA_full", "peak_SDA_max", "time_peak_SDA_max", "peak_SDA_mean", "time_peak_SDA_mean", "smr_type" )
  SDAdata_stepIntegral<-matrix(ncol=15, nrow=0)
  colnames(SDAdata_stepIntegral)<-c("ID","SMR","spar", "SDA_integrated", "end_SDA_estim_hr", "SMR_intergrated", "peak_SDA", "time_peak_SDA", "percentSMR_peak_SDA", "MO2_SDA_full", "peak_SDA_max", "time_peak_SDA_max", "peak_SDA_mean", "time_peak_SDA_mean", "smr_type" )

  # start of for loop applying the SDA function
  for(i in 1:length(unique(d_SMRsum$Ch))){

    if (length(unique(d_SMRsum$Ch))==1){
      Y0<-d_SMRsum
      Y0.full<-Ch.data.sda.full
    }else{
      Y0<-as.data.frame(Ch.data.sda[i])
      Y0.full<-as.data.frame(Ch.data.sda.full[i])
    }

    colnames(Y0)<-c("Ch",  "ID", "resp.V", "bw", "hour", "mo2_mean",  "mo2_sd", "mo2_min", "t_mean", "n", "mo2_max")
    colnames(Y0.full)<-c("Ch1.time_frame","min_start" ,"r2", "b" ,"m" ,"t_min", "Ch1.t_max","t_mean" ,"Ch" ,"DateTime_start", "Ch1.type","n_min", "Ch1.ID_code","bw","mo2","ID" ,"resp.V", "hour")

    mo2<-Y0$mo2_min
    BW<-Y0$bw[1]
    ID<-Y0$ID[1]
    N_mo2<-length(mo2)

    Y.Ch<-as.character(Y0$Ch[1])
    d<-Y0

    # peak SDA from different measurements
    # peak_SDA_max<-max(Y0.full$mo2)[1] # if more than one value is "max" then take the first one, that is also used for time calculations
    # time_peak_SDA_max<-Y0.full$min_start[which(Y0.full$mo2==peak_SDA_max)] # in minutes

    if(local_path | !dir.exists("SDA")){
      SDAplot_name <-	paste(data.SDA, "_", d$Ch[1], "_SDA_PLOT.png", sep='')
    }else{
      SDAplot_name <-	paste("./SDA/plots_ch_SDA/", data.SDA, "_", d$Ch[1], "_SDA_PLOT.png", sep='')
    }

    # Find the smr value to use for SDA calculations
    if (sda_threshold_level[1]=="SMR_vals" & SMR_calc==FALSE){
      # print(Y.Ch)
      b<-SMR_vals[as.numeric(substr(Y.Ch, start=3, stop=3))]
      b<-b*as.numeric(sda_threshold_level[2])

    }else{
      # either SMR_calc=TRUE in which case the newdata.smr will be the newly analyzed dataframe
      # or SMR_calc=FALSE in which case the newdata.smr will be the imported data frame (data.MR)
      smr.row<-newdata.smr[which(as.character(newdata.smr$Ch)==as.character(d$Ch[1])),]
      smr.row[,c(4:ncol(smr.row))] <- lapply(smr.row[,c(4:ncol(smr.row))], as.character)
      smr.row[,c(4:ncol(smr.row))] <- lapply(smr.row[,c(4:ncol(smr.row))], as.numeric)

      ID<-smr.row["ID"]

      # print(smr.row)
      if(sda_threshold_level[1]=="SMR_mean10minVal"){
        b<-as.numeric(round(smr.row["smr_mean10minVal"],2))
      }
      if(sda_threshold_level[1]=="SMR_low10quant"){
        b<-as.numeric(round(smr.row["SMR_low10quant"],2))
      }
      if(sda_threshold_level[1]=="SMR_low15quant"){
        b<-as.numeric(round(smr.row["SMR_low15quant"],2))
      }
      if(sda_threshold_level[1]=="SMR_low20quant"){
        b<-as.numeric(round(smr.row["SMR_low20quant"],2))
      }
      if(MLND & sda_threshold_level[1]=="smr_mlnd"){
        b<-as.numeric(round(smr.row["smr_mlnd"],2))
      }

      b<-b*as.numeric(sda_threshold_level[2])

    }

    end_SDA<-end_SDA_Ch[as.numeric(substr(Y.Ch, start=3, stop=3))]
    # messages and function checks
    if(!is.numeric(b) | !exists("b")){
      stop("provide usable SMR value to use for SDA calculations")
    }

    if(i==1){
      message(paste("SDA data: ",Y0$n[1], " MO2 measurements each hour;  ", N_mo2, " total hours", sep=""))
    }
    message(paste(d$Ch[1],": SMR value = ", b, sep =""))

    # Calculate the area under the curve using a point by point, direct integrals.
    # use the defined end SDA end value, or when it first hits the SMR defined value

    spars <- c(0.1,0.2,0.3)
    zero.row<-d[1,] # the first row of the data (experimentally typical)
    # print(Y.Ch)

    # d$hour<-d$hour+0.5 # make a half hour because the first measurement and the first hour mean are both 0, fitting smooth spline one will be dropped

    if(begin_smr_hr_zero){ # starting measurement at the beginning of the file
      # zero.row<-d[1,]
      zero.row$mo2_mean<-b # replacing with the selected SMR value
      zero.row$mo2_min<-b # replacing with the selected SMR value
      zero.row$hour<-0 # replacing with the selected SMR value

      d<-rbind(zero.row, d)

      if(i==1){
        d_SMRsum_wDELAY<-d
      }else{
        d_SMRsum_wDELAY<-rbind(d_SMRsum_wDELAY,d)
      }
    }

    if(nrow(d)>=4){
      for (n in 1:length(spars)){
        # if (b<1) {next}
        if (n == 1) {
					png(SDAplot_name, width = 4, height = length(spars)*3, units="in", res=200)
					par(mfrow=c(length(spars),1), mar=c(4,4,3,1)+0.5)
        }
        # spar,d, SDAdata, b, sda_threshold, end_SDA, begin_smr_hr_zero
        SDAdata <- SDA.spar(spar = spars[n], d, SDAdata = SDAdata, b, sda_threshold = sda_threshold_level[2], end_SDA = end_SDA, begin_smr_hr_zero = begin_smr_hr_zero)
        # ncol(SDAdata)
        if (n == length(spars)){
          dev.off()
        }
      }

    }else{
      message("Not enough points to do smoothing & integration for SDA (n < 4)")
    }

    if(is.na(end_SDA_Ch[i])){
      end_SDA_absolute<-d$hour[which(d$mo2_min[1:nrow(d)] <= b)[2]] # excludes the first hour since that is manually added to be SMR level
      end_row<-which(d$mo2_min[1:nrow(d)] <= b)[2]
      # 	   print(c(end_SDA_absolute, b, end_row))
      #      print(d[which(d$mo2_min[1:nrow(d)] <= b),])
      #
      if(is.na(end_row)){
        end_row<-nrow(d)
      }

    }else{
      end_SDA_absolute<-end_SDA_Ch[i]
      end_row<-which(d$hour== end_SDA_absolute)
    }
    message(paste("   SDA ends at ", d$hour[end_row], " h", sep=""))

    ### INTEGRATION CODE
    Full_SDA<-AUC(x=d$hour[1:end_row], y=d$mo2_min[1:end_row], method="trapezoid")
    SMRchunk<-AUC(x=d$hour[1:end_row], y=rep(b, end_row), method="trapezoid")

    peak_SDA<-max(d$mo2_min[1:end_row])[1]
    peak_SDA_max<-max(d$mo2_max[1:end_row])[1]
    peak_SDA_mean<-max(d$mo2_mean[1:end_row])[1]
    # }

    time_peak_SDA<-d$hour[which(d$mo2_min==peak_SDA)]
    time_peak_SDA_max<-d$hour[which(d$mo2_max==peak_SDA_max)[1]]
    time_peak_SDA_mean<-d$hour[which(d$mo2_mean==peak_SDA_mean)[1]] # in h not saved in the dataframe
    percentSMR_peak_SDA<-round((peak_SDA/b)*100,2)

    SDA_integral<- Full_SDA-SMRchunk

    integral_values<-as.data.frame(t(c(as.character(SDAdata$ID[nrow(SDAdata)]), b, "AUC",
                                       SDA_integral, end_SDA_absolute,SMRchunk,
                                       as.character(peak_SDA),
                                       as.character(time_peak_SDA),
                                       as.character(percentSMR_peak_SDA), NA,
                                       as.character(peak_SDA_max),
                                       as.character(time_peak_SDA_max),
                                       as.character(peak_SDA_mean),
                                       as.character(time_peak_SDA_mean),
                                       as.character(SDAdata$smr_type[nrow(SDAdata)]))))
    colnames(integral_values)<-colnames(SDAdata_stepIntegral)
    SDAdata_stepIntegral<-rbind(SDAdata_stepIntegral, integral_values)

  }# end of for loop applying the SDA function

  if(begin_smr_hr_zero){
    d_SMRsum_wDELAY<-merge(d_SMRsum_wDELAY, unique(SDAdata[,1:2]),  by.x="ID")
    d_SMRsum_wDELAY$SMR<-as.numeric(as.character(d_SMRsum_wDELAY$SMR))

    lst <- lapply(d_SMRsum_wDELAY , unlist)
    d_SMRsum_wDELAY <- (data.frame(lapply(lst, `length<-`, max(lengths(lst)))))
    write.csv(file=SDAhrlydata_name_wDELAY, d_SMRsum_wDELAY, row.names=FALSE)
  }

  d_SMRsum<-merge(d_SMRsum, unique(SDAdata[,1:2]),  by.x="ID")
  d_SMRsum$SMR<-as.numeric(as.character(d_SMRsum$SMR))

  lst <- lapply(d_SMRsum, unlist)
  d_SMRsum <- (data.frame(lapply(lst, `length<-`, max(lengths(lst)))))
  write.csv(file=SDAhrlydata_name, d_SMRsum, row.names=FALSE)

  d_SMRsum$ID<-as.character(d_SMRsum$ID)

  sda_hr_plot<-ggplot(data=d_SMRsum, aes(y=mo2_mean, x=hour))+
    geom_point(data=d_SMR_original, mapping = aes(y = mo2, x = min_start/60), pch=19, color = "grey", alpha = 0.5)+
    geom_point(size=2, pch=21, fill="grey30", alpha=0.9, colour="black")+
    geom_point(data=d_SMRsum, aes(y=mo2_min, x=hour), pch="-", size=2, stroke =3,color="black", alpha=0.7)+
    geom_line(data=d_SMRsum, aes(y=mo2_min, x=hour), size=1, alpha=0.7)+
    geom_errorbar(ymin=d_SMRsum$mo2_mean-d_SMRsum$mo2_sd, ymax = d_SMRsum$mo2_mean+d_SMRsum$mo2_sd, alpha=0.5 )+
    theme_classic()+
    geom_hline(aes(yintercept=SMR), data=d_SMRsum, lty=1)+
    # geom_hline(aes(yintercept=SMR*0.9), data=d_SMRsum, lty=2, colour="grey")+
    # geom_hline(aes(yintercept=SMR*1.1), data=d_SMRsum, lty=2, colour="grey")+
    ggtitle(paste(sda_threshold_level[1], " at (%): ", as.numeric(sda_threshold_level[2])*100, sep=""))+
    ylab("grey = MO2 mean +/- SE, black = hourly minimum recorded MO2")+
    # geom_points(aes(x=
    facet_wrap(.~ID, ncol=1, nrow=n_id, scales="free")
  if(begin_smr_hr_zero){
    # sda_hr_plot<- sda_hr_plot + geom_point(data=d_SMRsum_wDELAY, aes(y=mo2_min, x=hour), pch=21, size=1, fill="red", alpha=0.7)
    sda_hr_plot<- sda_hr_plot + geom_line(data=d_SMRsum_wDELAY, aes(y=mo2_min, x=hour), colour="red", alpha=0.7)
    # sda_hr_plot<- sda_hr_plot +	geom_point(data=plotDF, aes(y=mo2_min, x=hour), colour="green", pch=8)
  }

  png(plotname.sda.data, width=6, height=10, units="in", res=300)
    print(sda_hr_plot)
  dev.off()

  lst <- lapply(SDAdata, unlist)
  colnames(SDAdata_stepIntegral)<-colnames(SDAdata)
  SDAdata <- (data.frame(lapply(lst, `length<-`, max(lengths(lst)))))
  SDAdata<-rbind(SDAdata_stepIntegral, SDAdata)

  write.csv(file = SDAdata_name, SDAdata, row.names=FALSE)

}




