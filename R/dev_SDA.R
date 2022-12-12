# Nov 30 notes:
# add -- >
# channel specific delays (the same concept as MMR putting one fish in the channel at the time)
# specify channel specific start and stop times for the SDA portion.
# data.SMR= the raw csv file to be analyzed for SMR

SDA<-function(AnimalID,
              BW.animal,
              resp.V,
              r2_threshold_smr,
              sda_threshold_level,
              scaling_exponent_mmr = 1,
              scaling_exponent_smr = 1,
              date_format = c("%m/%d/%Y %H:%M:%S", "GMT"), # herehere
              data.SDA = NULL,
              analyzed_MR = NULL,
              SMR_calc = TRUE,
              SMR_vals = c(NA, NA, NA, NA),
              drop_ch = NULL,
              N_Ch = 4,
              end_SDA_Ch = NA,
              MLND = TRUE,
              background_prior = NULL, # herehere
              background_post = NULL,# herehere
              background_slope = NULL,
              background.V = NULL,
              match_background_Ch = FALSE,
              background_linear_gr = FALSE,
              local_path = FALSE, # herehere
              handling_delay = 0,
              begin_smr_hr_zero=FALSE){

  SDA.spar<-function(spar,d, SDAdata, b, sda_threshold, end_SDA, begin_smr_hr_zero){
  d$hour<-as.numeric(as.character(d$hour))

  fit<-smooth.spline(d$hour,d$mo2_min, spar=spar)
  f = function(x) {predict(fit, x)$y}

  end<-round(d$hour[nrow(d)],1)
  newx<-seq(0,end, by=0.1)
  newy<-predict(fit, newx, deriv=0)
  lapply(newy, as.numeric)
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


    if(is.na(end_SDA) & spar==0.1){
      if (begin_smr_hr_zero==TRUE){
        end_SDA <- newVal$init[(which(round(newVal$V2[2:nrow(newVal)], 3)<=b))[1]] # force the function to look beyond te first value, which is set to SMR level
      }else{
        end_SDA<-newVal$init[nrow(newVal)]
        message("Respiration post digestion does not reach the chosen SMR levels: ",b, "  - suggest changing the sda_threshold")
      }

    }

    if(is.na(end_SDA)){
      end_SDA<-newVal$init[nrow(newVal)]
    }
  }


  if(end_SDA > max(newVal$init)){
    message("The SDA end value is greater than the duration of respirometry trial | end_SDA set to the last value")
    end_SDA<-newVal$init[nrow(newVal)]
  }

  #3 integrate recovery curve with to the provided end EPOC time
  f.int = function(x) {integrate(f, lower=0, upper=end_SDA)$value}
  f.int.vec = Vectorize(f.int, vectorize.args='x')
  f.int.vec = Vectorize(f.int, vectorize.args='x')
  full<-f.int.vec(end_SDA)
  # SMR block
  SMR<-integrate(f.smr, lower=0, upper=end_SDA)$value
  SDA_full<-round(full-SMR,3) # the costs of digestion
  MO2_full<-round(newVal$V2[which(newVal$init==end_SDA)],3)

  end_SDA_row<-which(d$hour==round(end_SDA,0))
  # print(c(end_SDA_row, end_SDA))
  # print(d[end_SDA_row,])
  # find peak SDA - use min values

  # 		if(end_SDA == d$hour[nrow(d)]){
  #   		peak_SDA<-max(d$mo2_min)
  #   		peak_SDA_max<-max(d$mo2_max)
  #   		peak_SDA_mean<-max(d$mo2_mean)
  # }else{
  peak_SDA<-max(d$mo2_min[1:end_SDA_row])[1]
  # print(max(d$mo2_min[1:end_SDA_row])[1])
  peak_SDA_max<-max(d$mo2_max[1:end_SDA_row])[1]
  peak_SDA_mean<-max(d$mo2_mean[1:end_SDA_row])[1]
  # }

  time_peak_SDA<-d$hour[which(d$mo2_min==peak_SDA)]
  time_peak_SDA_max<-d$hour[which(d$mo2_max==peak_SDA_max)[1]]
  time_peak_SDA_mean<-d$hour[which(d$mo2_mean==peak_SDA_mean)[1]] # in h not saved in the dataframe
  percentSMR_peak_SDA<-round((peak_SDA/b)*100,2)

  values<-as.data.frame(t(c(as.character(d$ID[1]),
                            b, spar, SDA_full, end_SDA,
                            SMR, peak_SDA, time_peak_SDA, percentSMR_peak_SDA, MO2_full,
                            peak_SDA_max, time_peak_SDA_max, peak_SDA_mean, time_peak_SDA_mean,  sda_threshold)))



  # print(c(peak_SDA_max, peak_SDA_mean, peak_SDA, time_peak_SDA_mean, time_peak_SDA, time_peak_SDA_max, SMR, spar))
  # print(c(b, MO2_full, SDA_full, sda_threshold, percentSMR_peak_SDA))
  # if MO2_full and b are the same(ish - that's indicating that something is not right)
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
  legend("topright", legend=paste("spar = ", spar))
  # text(x=(max(d$mo2_min)+50)-(0.01*(max(d$hour)+50)),y=(max(d$mo2_min)+2)-(scale[i]*(max(d$mo2_min)+2)), label=paste("EPOC=",EPOC_full,"/ ",smr_type, sep=""), cex=0.8, col=col_smr[i], pos=2)
  points(x=time_peak_SDA, y=peak_SDA, pch=23, col="red",cex=2)
  points(x=time_peak_SDA_mean, y=peak_SDA_mean, pch=23, col="blue",cex=2)
  # text(x=(max(d$time_mo2)+50)-(0.01*(max(d$time_mo2)+50)),y=(max(d$mo2)+2)-(scale[i]*(max(d$mo2)+2)), label=paste("EPOC=",EPOC_full,"/ ",smr_type, sep=""), cex=0.8, col=col_smr[i], pos=2)


  SDAdata<-rbind(SDAdata,SDAdata.temp)

  return(SDAdata)

} # end SDA.spar

  graphics.off()

  filename.SMR<-paste(gsub('.{4}$', '',data.SDA), "_SMR", sep="")

  # START -- >>> SMR data to be calculated from the file
  # if (SMR_calc==TRUE){
  # START-- >>> background
  # ----##


  # background file manipulation
  # 1. find what channels recorded background
  if (!is.na(background_post) | !is.na(background_prior) ) {
    # 2 calculate mean for each background channel
    if((is.null(background.V) & !is.null(background_slope) ) |(!is.null(background.V) & is.null(background_slope))) {
      stop_function <- TRUE
      if(stop_function) stop("If using manually input background slopes, must provide both, a unique slope and volume of the repirometer")
    }


    back_prior<-read.csv(background_prior)
    back_post<-read.csv(background_post)
    back_ch<-length (unique(back_post$Ch))

    # Jan 4 2020: make linear regression over time and correct background based on a predicted value
    # create
    if (background_linear_gr==TRUE){
      if(!exists("back_prior") | !exists("back_post")){
        stop("Missing a file to model the growth of bacteria, must provide both: \"back_prior\" \"and back_post\" files")
      }

      if (!dir.exists("../plots_background")){
        dir.create(file.path("../plots_background"), recursive = TRUE)
      }

      # message("Background: assuming a linear growth of bacteria over time | Using channel specific growth slopes")

      # getting the right date and time format
      back_prior$DateTime_start<-as.character(gsub("(", "", back_prior$DateTime_start, fixed=TRUE))
      back_prior$DateTime_start<-as.character(gsub(")", "", back_prior$DateTime_start, fixed=TRUE))
      back_post$DateTime_start<-as.character(gsub("(", "", back_post$DateTime_start, fixed=TRUE))
      back_post$DateTime_start<-as.character(gsub(")", "", back_post$DateTime_start, fixed=TRUE))

      if (date_format== "m/d/y"){
        back_prior$DateTime_start<- strptime(back_prior$DateTime_start, format="%m/%d/%y %H:%M:%OS")
        back_post$DateTime_start<- strptime(back_post$DateTime_start, format="%m/%d/%y %H:%M:%OS")
      }
      if (date_format== "d/m/y"){
        back_prior$DateTime_start<- strptime(back_prior$DateTime_start, format="%d/%m/%y %H:%M:%OS")
        back_post$DateTime_start<- strptime(back_post$DateTime_start, format="%d/%m/%y %H:%M:%OS")
      }
      if (date_format== "y-m-d"){
        # DateTime<- chron(dates=back_$date,times=back_$time,format=c('y-m-d','h:m:s'))
        back_prior$DateTime_start<- strptime(back_prior$DateTime_start, format="%y-%m-%d %H:%M:%OS")
        back_post$DateTime_start<- strptime(back_post$DateTime_start, format="%y-%m-%d %H:%M:%OS")
      }

      back_all<- rbind(back_prior, back_post)
      back_all$DateTime_start<- as.POSIXct(back_all$DateTime_start)

      back_regression_plot<-ggplot(data=back_all, aes(x=DateTime_start, y=m, colour=Ch, group=Ch))+
        geom_point()+
        geom_smooth(method="lm", se=FALSE)+
        theme_bw()+
        theme(axis.text.x = element_text(angle = 45))+
        facet_grid(Ch~.)

      if (path == "."){
        plotname.backgr.analysis<-paste( filename.SMR,"_PLOT_BACKGROUND_regressions.png", sep="")
        plotname.backgr<-paste( filename.SMR,"_PLOT_BACKGROUND.png", sep="")
      }else{
        plotname.backgr.analysis<-paste("../plots_background/", filename.SMR,"_PLOT_BACKGROUND_regressions.png", sep="")
        plotname.backgr<-paste("../plots_background/", filename.SMR,"_PLOT_BACKGROUND.png", sep="")

      }
      png(plotname.backgr.analysis, width=4, height=7, res=300, units="in")
      print(back_regression_plot)
      dev.off()


      if (match_background_Ch==TRUE){
        back_ch_regressions<-list()

        for(i in 1:back_ch){
          back_ch_d<- back_all[back_all$Ch==(unique(back_all$Ch))[i],]
          Ch<-substr(as.character(back_ch_d$Ch[1]), start=3, stop=3)

          back_regression_name<- paste("back_regression", Ch, sep="") # channel names with a channel # at the end
          regression<- lm(m~DateTime_start, data = back_ch_d)
          assign(back_regression_name, regression)
          back_ch_regressions[[i]] <- assign(back_regression_name, regression)

        }
        # WORK
        # save background regressions

      }else{# end for ch specific regressions
        # get one background slope based on all datapoints collected
        back_regression<- lm(m~DateTime_start, data = back_all)
        message("Background: assuming a linear growth of bacteria over time | Using one mean slope for all channels")
      }


    } # end for getting linear regressions for the background

    if (background_linear_gr==FALSE){
      if (!is.na(background_prior)){
        back_ch_prior<-list()
        back_ch_prior_names<-list()

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
      if (!is.na(background_post)){

        back_ch_post<-list()
        back_ch_post_names<-list()

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

      ## should look for possible variables:
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


  newdata.smr<-as.data.frame(matrix(ncol=17, nrow=0))
  names(newdata.smr)<-c("filename", "ID", "Ch", "BW","t_min","t_max", "t_mean", "N_mo2", #8
                        "smr_mean10minVal","smr_SD10minVal", "smr_CV10minVal", "SMR_low10quant","SMR_low15quant","SMR_low20quant", #6
                        "smr_mlnd", "smr_CVmlnd", "smr_Nmlnd")#3

  cols = c(4:17)
  newdata.smr[,cols] %<>% lapply(function(x) as.numeric(as.character(x)))
  cols2 = c(1:3)
  newdata.smr[,cols2] %<>% lapply(function(x) as.character(x))

  # print(ncol(newdata.smr))

  d_SMR<-read.csv(data.SDA)
  # drop any unwanted channels
  if(!is.null(drop_ch[1])){
    n_ch_drop<-length(drop_ch)
    for(i in 1:n_ch_drop){
      d_SMR<-d_SMR[(substr(as.character(d_SMR$Ch), start=3, stop=3)!=drop_ch[i]),]
    }
  }
  d_SMR$Ch<-factor(d_SMR$Ch)


  if (path == "."){
    plotname.freq<-paste( filename.SMR,"_PLOT_SMR_analyses.png", sep="")
    plotname.smr.meth<-paste( filename.SMR,"_PLOT_SMR_methodsALL.png", sep="")

  }else{
    plotname.freq<-paste("../plots_min_values_SMR/", filename.SMR,"_PLOT_SMR_analyses.png", sep="")
    plotname.smr.meth<-paste("../plots_methods_sum_SMR/", filename.SMR,"_PLOT_SMR_methodsALL.png", sep="")
  }


  if(!colnames(d_SMR)[11]=="type"){
    d_SMR$type="SMR"
    d_SMR<-d_SMR[,c("time_frame", "min_start", "r2", "b", "m", "t_min", "t_max", "t_mean" ,"Ch", "DateTime_start", "type", "n_min", "ID_code" )]
  }


  ## choose what to keep and what not for the SMR
  # 2) keep only > 2 min sections for SMR calculations any type  selected above (SMR, pre-, post-shallow slopes) and exclude "SMR-cut", "SMR-cut1", "SMR-cut2"
  d_SMR.1<-d_SMR[d_SMR$n_min>=1,]
  if (nrow(d_SMR.1)==0){
    warning("!! SMR measurements shorter than 1 min !!")
  }else{
    d_SMR<-d_SMR[d_SMR$n_min>=1,]
    message("Not using SMR measurements < 1 min")
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


  d_SMR$DateTime_start<-as.character(gsub("(", "", d_SMR$DateTime_start, fixed=TRUE))
  d_SMR$DateTime_start<-as.character(gsub(")", "", d_SMR$DateTime_start, fixed=TRUE))

  if (date_format== "m/d/y"){
    d_SMR$DateTime_start<- strptime(d_SMR$DateTime_start, format="%m/%d/%y %H:%M:%OS")
  }
  if (date_format== "d/m/y"){
    d_SMR$DateTime_start<- strptime(d_SMR$DateTime_start, format="%d/%m/%y %H:%M:%OS")
  }
  if (date_format== "y-m-d"){
    # DateTime<- chron(dates=back_$date,times=back_$time,format=c('y-m-d','h:m:s'))
    d_SMR$DateTime_start<- strptime(d_SMR$DateTime_start, format="%y-%m-%d %H:%M:%OS")
  }


  # 000 Jan 4 addition - account for background using linear regression
  # making predictions
  if(background_linear_gr==TRUE & match_background_Ch==FALSE){
    background_slopes<-data.frame(DateTime_start = d_SMR$DateTime_start)
    background_slopes$back_m<-predict(back_regression, background_slopes)

    if (path == "."){
      plotname.backgr.analysis<-paste( filename.SMR,"_PLOT_BACKGROUND_regressions.png", sep="")
      plotname.backgr<-paste( filename.SMR,"_PLOT_BACKGROUND.png", sep="")
    }else{
      plotname.backgr.analysis<-paste("../plots_background/", filename.SMR,"_PLOT_BACKGROUND_regressions.png", sep="")
      plotname.backgr<-paste("../plots_background/", filename.SMR,"_PLOT_BACKGROUND.png", sep="")

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

    d_SMR$background_slope<- background_slopes$back_m

  }



  if(background_linear_gr==TRUE & match_background_Ch==TRUE){

    message("SMR corrected for background: using Ch specific average background")
    background_slopes<-data.frame(DateTime_start = d_SMR$DateTime_start)

    for (i in 1:nrow(d_SMR)){

      # print(c(d_SMR$bw[1],d_SMR$Ch[i]))
      if(substr(d_SMR$Ch[i], start=3, stop=3) == "1"){
        background_slopes$back_m[i]<-predict(back_regression1, data.frame(DateTime_start = d_SMR$DateTime_start[i]))
        back_m1<-predict(back_regression1, data.frame(DateTime_start = d_SMR$DateTime_start[i]))
        d_SMR$mo2[i]<-(( d_SMR$m[i]*(d_SMR$resp.V[i]-d_SMR$bw[i])) - (back_m1 * d_SMR$resp.V[i])) /(d_SMR$bw[i])#^scaling_exponent_smr) # units mgO2 kg-1 min-1 -CORRECTED for back resp
        # message("Channel: 1,  bacterial respiration slope: ", back_m1)
      }

      if(substr(d_SMR$Ch[i], start=3, stop=3) == "2"){
        background_slopes$back_m[i]<-predict(back_regression2, data.frame(DateTime_start = d_SMR$DateTime_start[i]))
        back_m2<-predict(back_regression2, data.frame(DateTime_start = d_SMR$DateTime_start[i]))
        d_SMR$mo2[i]<-(( d_SMR$m[i]*(d_SMR$resp.V[i]-d_SMR$bw[i])) - (back_m2 * d_SMR$resp.V[i])) /(d_SMR$bw[i])#^scaling_exponent_smr) # units mgO2 kg-1 min-1 -CORRECTED for back resp
        # message("Channel: 2,  bacterial respiration slope: ", back_m2)
      }

      if(substr(d_SMR$Ch[i], start=3, stop=3) == "3"){
        background_slopes$back_m[i]<-predict(back_regression3, data.frame(DateTime_start = d_SMR$DateTime_start[i]))
        back_m3<-predict(back_regression3, data.frame(DateTime_start = d_SMR$DateTime_start[i]))
        d_SMR$mo2[i]<-(( d_SMR$m[i]*(d_SMR$resp.V[i]-d_SMR$bw[i])) - (back_m3 * d_SMR$resp.V[i])) /(d_SMR$bw[i])#^scaling_exponent_smr) # units mgO2 kg-1 min-1 -CORRECTED for back resp
        # message("Channel: 3,  bacterial respiration slope: ", back_m3)
      }

      if(substr(d_SMR$Ch[i], start=3, stop=3) == "4"){
        background_slopes$back_m[i]<-predict(back_regression4, data.frame(DateTime_start = d_SMR$DateTime_start[i]))
        back_m4<-predict(back_regression4, data.frame(DateTime_start = d_SMR$DateTime_start[i]))
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
  if ((( !is.na(background_post) | !is.na(background_prior)) & match_background_Ch==FALSE) & is.null(background_slope) & background_linear_gr==FALSE){
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
  if ((!is.na(background_post) | !is.na(background_prior)) & match_background_Ch==TRUE & background_linear_gr==FALSE){
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
  if ((is.na(background_post) & is.na(background_prior)) & is.null(background_slope)){
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
    ggtitle(data.SDA)

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
    geom_line(size=0.5, alpha=0.)+
    theme_classic()+
    ggtitle("lowest 10 VALUES  excluding 5 smallest")+
    theme(legend.position="top")+
    facet_grid(Ch~.)

  # the 10% percentile
  # perc=c(0.1,0.15, 0.2)

  a2<-as.data.frame(matrix(ncol=3, nrow=0))

  for (i in unique(d_SMR$Ch)){
    split_temp<-as.data.frame(split(d_SMR, d_SMR$Ch)[i])
    colnames(split_temp)<-colnames(d_SMR)
    quant10<-quantile(split_temp$mo2, 0.1, na.rm=FALSE)
    quant15<-quantile(split_temp$mo2, 0.15, na.rm=FALSE)
    quant20<-quantile(split_temp$mo2, 0.2, na.rm=FALSE)

    # colnames(a2)<-c("Ch","quantiles", "mo2_perc")
    a2<-rbind(a2,as.data.frame(t(c(as.character(split_temp$Ch[1]), "10%", as.numeric(quant10)))))
    # colnames(a2)<-c("Ch","quantiles", "mo2_perc")
    a2<-rbind(a2,as.data.frame(t(c(as.character(split_temp$Ch[1]), "15%", as.numeric(quant15)))))
    # colnames(a2)<-c("Ch","quantiles", "mo2_perc")
    a2<-rbind(a2,as.data.frame(t(c(as.character(split_temp$Ch[1]), "20%", as.numeric(quant20)))))
  }

  colnames(a2)<-c("Ch","quantiles", "mo2_perc")
  a2$mo2_perc<-as.numeric(as.character( a2$mo2_perc))

  #
  # 		a2<-d_SMR[,2:ncol(d_SMR)] %>%
  # 			group_by(Ch)%>%
  # 			summarise( quantiles = list(sprintf("%1.0f%%", perc*100)),
  # 			mo2_perc = list(quantile(mo2, perc, na.rm=FALSE))) %>%
  # 			unnest(cols = c(quantiles, mo2_perc))
  #
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

  min10_percPlot<-ggplot(data=d_SMR, aes(x=min_start, y=mo2))+
    geom_point(size=1)+
    geom_point(data=d_SMR[which(d_SMR$Ch=="Ch1" & d_SMR$mo2<=a2$mo2_perc[row_ch1_10perc]),], aes(x=min_start, y=mo2), colour="blue",size=3)+
    geom_point(data=d_SMR[which(d_SMR$Ch=="Ch2" & d_SMR$mo2<=a2$mo2_perc[row_ch2_10perc]),], aes(x=min_start, y=mo2), colour="blue",size=3)+
    geom_point(data=d_SMR[which(d_SMR$Ch=="Ch3" & d_SMR$mo2<=a2$mo2_perc[row_ch3_10perc]),], aes(x=min_start, y=mo2), colour="blue",size=3)+
    geom_point(data=d_SMR[which(d_SMR$Ch=="Ch4" & d_SMR$mo2<=a2$mo2_perc[row_ch4_10perc]),], aes(x=min_start, y=mo2), colour="blue",size=3)+
    geom_line(size=0.5, alpha=0.5)+
    theme_classic()+
    ggtitle("10 PERCENTILE")+
    theme(legend.position="top")+
    facet_grid(Ch~.)

  min15_percPlot<-ggplot(data=d_SMR, aes(x=min_start, y=mo2))+
    geom_point(size=1)+
    geom_point(data=d_SMR[which(d_SMR$Ch=="Ch1" & d_SMR$mo2<=a2$mo2_perc[row_ch1_15perc]),], aes(x=min_start, y=mo2), colour="darkturquoise",size=3)+
    geom_point(data=d_SMR[which(d_SMR$Ch=="Ch2" & d_SMR$mo2<=a2$mo2_perc[row_ch2_15perc]),], aes(x=min_start, y=mo2), colour="darkturquoise",size=3)+
    geom_point(data=d_SMR[which(d_SMR$Ch=="Ch3" & d_SMR$mo2<=a2$mo2_perc[row_ch3_15perc]),], aes(x=min_start, y=mo2), colour="darkturquoise",size=3)+
    geom_point(data=d_SMR[which(d_SMR$Ch=="Ch4" & d_SMR$mo2<=a2$mo2_perc[row_ch4_15perc]),], aes(x=min_start, y=mo2), colour="darkturquoise",size=3)+
    geom_line(size=0.5, alpha=0.5)+
    theme_classic()+
    ggtitle("15 PERCENTILE")+
    theme(legend.position="top")+
    facet_grid(Ch~.)

  min20_percPlot<-ggplot(data=d_SMR, aes(x=min_start, y=mo2))+
    geom_point(size=1)+
    geom_point(data=d_SMR[which(d_SMR$Ch=="Ch1" & d_SMR$mo2<=a2$mo2_perc[row_ch1_20perc]),], aes(x=min_start, y=mo2), colour="deepskyblue2",size=3)+
    geom_point(data=d_SMR[which(d_SMR$Ch=="Ch2" & d_SMR$mo2<=a2$mo2_perc[row_ch2_20perc]),], aes(x=min_start, y=mo2), colour="deepskyblue2",size=3)+
    geom_point(data=d_SMR[which(d_SMR$Ch=="Ch3" & d_SMR$mo2<=a2$mo2_perc[row_ch3_20perc]),], aes(x=min_start, y=mo2), colour="deepskyblue2",size=3)+
    geom_point(data=d_SMR[which(d_SMR$Ch=="Ch4" & d_SMR$mo2<=a2$mo2_perc[row_ch4_20perc]),], aes(x=min_start, y=mo2), colour="deepskyblue2",size=3)+
    geom_line(size=0.5, alpha=0.5)+
    theme_classic()+
    ggtitle("20 PERCENTILE")+
    theme(legend.position="top")+
    facet_grid(Ch~.)

  png(plotname.freq, width = 15, height = 10, units="in", res=200)
  grid.arrange(p_freq, min10_plot, min10_percPlot, min15_percPlot, ncol=2, nrow=2)
  dev.off()

  # MLND algorhytm from suppl JFB issue 88 Chabot, Steffensen, and Farrell 2016
  # SMR dataframes split based on the channel
  if (length(unique(d_SMR$Ch))>1){
    Ch.data.smr<-split(d_SMR, d_SMR$Ch)
  }else{
    Ch.data.smr<-1
  }


  for(i in 1:length(unique(d_SMR$Ch))){

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
    if (path == "."){
      plotname.mlnd<-paste( filename.SMR,"_", Y.Ch, "_SMR_analyzed.png", sep="")

    }else{
      plotname.mlnd<-paste("../plots_mlnd_SMR/", filename.SMR,"_", Y.Ch, "_SMR_analyzed.png", sep="")
    }


    # START -->>> MLND
    # loop for MLND function
    if(MLND==TRUE){

      mlnd.m1<-Mclust(mo2, G=1:4)
      m1<-densityMclust(mo2, verbose=TRUE) # same as mlnd
      clas.m1<-m1$classification #
      clas.m1.t <- as.data.frame(table(clas.m1))
      clas.m1.t$clas.m1 <- as.numeric(levels(clas.m1.t$clas.m1))
      boot1 <- MclustBootstrap(m1, nboot = 1000, type = "bs")


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
  if (path == "."){
    filename.smr<-paste( gsub('.{12}$', '', data.SDA), "SMR_analyzed.csv", sep='')
    filename.MR<-paste( gsub('.{12}$', '', data.SDA), "MR_analyzed.csv", sep='')
  }else{
    filename.smr<-paste("../csv_analyzed_SMR/",gsub('.{12}$', '', data.SDA), "SMR_analyzed.csv", sep='')
    filename.MR<-paste("../csv_analyzed_MR/", gsub('.{12}$', '', data.SDA), "MR_analyzed.csv", sep='')
  }


  lst <- lapply(newdata.smr, unlist)
  newdata.smr <- (data.frame(lapply(lst, `length<-`, max(lengths(lst)))))

  write.csv(file=filename.smr, d_SMR, row.names=FALSE)
  write.csv(file=filename.MR, newdata.smr, row.names=FALSE)
  message("Save SMR, and MR files")
  # print(head(newdata.smr))
  # print(head(d_SMR))

  # return(newdata.smr)

  # --- >>> end of SMR_clac == TRUE


  if (SMR_calc==FALSE & (is.null(analyzed_MR))){
    stop("Need to provide analyzed MR file or specific SMR values")
  }


  if (SMR_calc==FALSE){
    newdata.smr<-read.csv(analyzed_MR)
    message("Use SMR values from previosly analyzed data")
  }

  # SDA
  # 1. get the mean of each hour
  d_SMR$hour<-as.factor(floor(d_SMR$min_start/60))
  # d_SMR$hour[which(as.numeric(as.character(d_SMR$hour))<1)]<-0
  # print(str(d_SMR))
  d_SMR$ID<-as.factor(as.character(d_SMR$ID))
  # d_SMR$resp.V<-as.factor(as.character(d_SMR$resp.V))
  # d_SMR$bw<-as.factor(as.character(d_SMR$bw))
  # d_SMR$hour<-as.factor(as.character(d_SMR$hour))

  # print(d_SMR$hour)

  d_SMRsum<-d_SMR %>%
    group_by(Ch, ID, resp.V, bw, hour) %>%
    summarize(mo2_mean = mean(mo2, na.rm=TRUE), mo2_sd = sd(mo2, na.rm=TRUE), mo2_min = min(mo2, na.rm=TRUE), t_mean=mean(t_mean), n=length (mo2), mo2_max = max(mo2, na.rm=TRUE))

  # d_SMRsum<-d_SMR %>%
  #   group_by(Ch,hour) %>%
  #   select(Ch, hour, ID, resp.V, bw, mo2, t_mean) %>%
  #   summarize(mo2_mean = mean(mo2, na.rm=TRUE), mo2_sd = sd(mo2, na.rm=TRUE), mo2_min = min(mo2, na.rm=TRUE), t_mean=mean(t_mean), n=length (mo2), mo2_max = max(mo2, na.rm=TRUE))


  # print(d_SMRsum)

  d_SMRsum$hour<-as.numeric(as.character(d_SMRsum$hour))
  d_SMRsum$ID<-as.factor(as.character(d_SMRsum$ID))
  d_SMRsum$resp.V<-as.numeric(as.character(d_SMRsum$resp.V))
  d_SMRsum$bw<-as.numeric(as.character(d_SMRsum$bw))

  n_id<-length(unique(levels( d_SMRsum$ID)))
  # d_SMRsum$hour<-as.numeric(as.character(d_SMRsum$hour))


  if (path == "."){
    plotname.sda.data<-	paste(gsub('.{4}$', '', data.SDA),"_SDA_hourly_PLOT.png", sep='')
    SDAdata_name<-paste(gsub('.{4}$', '', data.SDA),"_SDA_analyzed.csv", sep='')
    SDAhrlydata_name<-paste(gsub('.{4}$', '', data.SDA),"_SDA_hrly_analyzed.csv", sep='')
    SDAhrlydata_name_wDELAY<-paste(gsub('.{4}$', '', data.SDA),"_SDA_hrly_wDELAY_analyzed.csv", sep='')


  }else{
    plotname.sda.data<-	paste("../plots_SDA_hourly/",gsub('.{4}$', '', data.SDA),"_SDA_hourly_PLOT.png", sep='')
    SDAdata_name<-paste("../csv_analyzed_SDA/", gsub('.{4}$', '', data.SDA),"_SDA_analyzed.csv", sep='')
    SDAhrlydata_name<-paste("../csv_analyzed_SDA_hrly/", gsub('.{4}$', '', data.SDA),"_SDA_hrly_analyzed.csv", sep='')
    SDAhrlydata_name_wDELAY<-paste("../csv_analyzed_SDA_hrly/", gsub('.{4}$', '', data.SDA),"_SDA_hrly_wDELAY_analyzed.csv", sep='')
  }


  # get split dataframes for data analysis of each individual
  if (length(unique(d_SMRsum$Ch))>1){
    Ch.data.sda<-split(d_SMRsum, d_SMRsum$Ch)
    Ch.data.sda.full<-split(d_SMR, d_SMR$Ch)
  }else{
    Ch.data.sda<-d_SMRsum
    Ch.data.sda.full<-d_SMR
  }

  # SDA analysis function on hourly data frame
  #######



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

    if (path == "."){
      SDAplot_name <-	paste(data.SDA, "_", d$Ch[1], "_SDA_PLOT.png", sep='')
    }else{
      SDAplot_name <-	paste("../plots_ch_SDA/", data.SDA, "_", d$Ch[1], "_SDA_PLOT.png", sep='')
    }

    # Find the smr value to use for SDA calculations
    if (sda_threshold_level[1]=="SMR_vals" & SMR_calc==FALSE){
      b<-SMR_vals[as.numeric(substr(Y.Ch, start=3, stop=3))]
      b<-b*as.numeric(sda_threshold_level[2])

    }else{
      # either SMR_calc=TRUE in which case the newdata.smr will be the newly analyzed dataframe
      # or SMR_calc=FALSE in which case the newdata.smr will be the imported data frame (data.MR)
      smr.row<-newdata.smr[which(as.character(newdata.smr$Ch)==as.character(d$Ch[1])),]
      smr.row[,c(4:ncol(smr.row))] %<>% lapply(function(x) as.numeric(as.character(x)))

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
      if(MLND == TRUE & sda_threshold_level[1]=="smr_mlnd"){
        b<-as.numeric(round(smr.row["smr_mlnd"],2))
      }

      b<-b*as.numeric(sda_threshold_level[2])

    }

    end_SDA<-end_SDA_Ch[as.numeric(substr(Y.Ch, start=3, stop=3))]
    # messages and function checks
    if(!is.numeric(b) | !exists("b")){
      stop("provide usable SMR value to use for SDA calculations")
    }

    if(i==1){message(paste("SDA data: ",Y0$n[1], " MO2 measurements each hour;  ", N_mo2, " total hours", sep=""))}
    message(paste("SMR value for ", d$Ch[1], " is: ", b))
    ##

    ### Calculate the area under the curve using a point by point, direct integrals.
    # use the defined end SDA end value, or when it first hits the SMR defined value

    spars <- seq(0.1,1, by=0.1)
    zero.row<-d[1,]
    d<-d[d$hour>=handling_delay,]
    if(begin_smr_hr_zero==TRUE){
      # zero.row<-d[1,]
      zero.row$mo2_mean<-b # replacing with the selected SMR value
      zero.row$mo2_min<-b # replacing with the selected SMR value
      zero.row$hour<-0 # replacing with the selected SMR value

      # print(head(d))

      d<-rbind(zero.row, d)
      if(i==1){
        d_SMRsum_wDELAY<-d
      }else{
        d_SMRsum_wDELAY<-rbind( d_SMRsum_wDELAY,d)
      }

    }

    if(nrow(d) == 4 || nrow(d)>4){
      for (n in 1:length(spars)){
        # if (b<1) {next}
        if (n == 1) {
          png(SDAplot_name, width = 8, height = 12, units="in", res=200)
          par(mfrow=c(5,2), mar=c(2,4,2,1)+0.1)
        }
        # print(c(spars[n], d, SDAdata, b, sda_threshold_level[1], end_SDA, begin_smr_hr_zero))
        SDAdata <- SDA.spar(spars[n], d, SDAdata, b, sda_threshold_level[1], end_SDA, begin_smr_hr_zero)
        ncol(SDAdata)
        if (n == length(spars)){
          # print(SDAdata)
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
    message(paste("SDA ends at hour: ", d$hour[end_row], sep=""))

    #
    # print(SDAdata)
    # print(str(SDAdata))
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


  if(begin_smr_hr_zero==TRUE){
    d_SMRsum_wDELAY<-merge(d_SMRsum_wDELAY, unique(SDAdata[,1:2]),  by.x="ID")
    d_SMRsum_wDELAY$SMR<-as.numeric(as.character(d_SMRsum_wDELAY$SMR))
  }

  d_SMRsum<-merge(d_SMRsum, unique(SDAdata[,1:2]),  by.x="ID")
  d_SMRsum$SMR<-as.numeric(as.character(d_SMRsum$SMR))

  lst <- lapply(d_SMRsum, unlist)
  d_SMRsum <- (data.frame(lapply(lst, `length<-`, max(lengths(lst)))))
  write.csv(file=SDAhrlydata_name, d_SMRsum, row.names=FALSE)

  lst <- lapply(d_SMRsum_wDELAY , unlist)
  d_SMRsum_wDELAY <- (data.frame(lapply(lst, `length<-`, max(lengths(lst)))))
  write.csv(file=SDAhrlydata_name_wDELAY, d_SMRsum_wDELAY, row.names=FALSE)

  # 	  print(SDAdata_stepIntegral)

  plotDF<-SDAdata_stepIntegral
  plotDF$mo2_min<-NA ### add mo2values
  # print(colnames(d_SMR))
  colnames(plotDF)[names(plotDF) == "end_SDA_estim_hr"]<-"hour"
  plotDF$hour<-as.numeric(as.character(plotDF$hour))

  plotDF$ID<-as.character(plotDF$ID)
  d_SMRsum$ID<-as.character(d_SMRsum$ID)

  plotDF<- plotDF %>% left_join(d_SMRsum[,c("ID", "mo2_min", "hour")], by= c("ID", "hour"))
  colnames(plotDF)  <- c("ID","SMR","spar", "SDA_integrated", "hour", "SMR_intergrated", "peak_SDA", "time_peak_SDA", "percentSMR_peak_SDA", "MO2_SDA_full", "peak_SDA_max", "time_peak_SDA_max", "peak_SDA_mean", "time_peak_SDA_mean", "smr_type", "" , "mo2_min")
  plotDF$hour<-as.numeric(as.character(plotDF$hour))
  plotDF$mo2_min<-as.numeric(as.character(plotDF$mo2_min))

  sda_hr_plot<-ggplot(data=d_SMRsum, aes(y=mo2_mean, x=hour))+
    geom_point(size=2, pch=21, fill="grey", alpha=0.9, colour="black")+
    geom_point(data=d_SMRsum, aes(y=mo2_min, x=hour), pch=21, size=3, fill="black", alpha=0.7)+
    geom_line(data=d_SMRsum, aes(y=mo2_min, x=hour), size=1, alpha=0.7)+
    geom_errorbar(ymin=d_SMRsum$mo2_mean-d_SMRsum$mo2_sd, ymax = d_SMRsum$mo2_mean+d_SMRsum$mo2_sd, alpha=0.5 )+
    theme_classic()+
    geom_hline(aes(yintercept=SMR), data=d_SMRsum, lty=1)+
    geom_hline(aes(yintercept=SMR*0.9), data=d_SMRsum, lty=2, colour="grey")+
    geom_hline(aes(yintercept=SMR*1.1), data=d_SMRsum, lty=2, colour="grey")+
    ggtitle(paste(sda_threshold_level[1], " proportionally adjusted to the level (%): ", as.numeric(sda_threshold_level[2])*100, sep=""))+
    ylab("grey = MO2 mean +/- SEM, black = hourly MO2 min ")+
    # geom_points(aes(x=
    facet_wrap(.~ID, ncol=1, nrow=n_id, scales="free")
  if(begin_smr_hr_zero==TRUE){
    sda_hr_plot<- sda_hr_plot + geom_point(data=d_SMRsum_wDELAY, aes(y=mo2_min, x=hour), pch=21, size=1, fill="red", alpha=0.7)
    sda_hr_plot<- sda_hr_plot + geom_line(data=d_SMRsum_wDELAY, aes(y=mo2_min, x=hour), colour="red", alpha=0.7)
    sda_hr_plot<- sda_hr_plot +	geom_point(data=plotDF, aes(y=mo2_min, x=hour), colour="green", pch=8)
  }

  png(plotname.sda.data, width=6, height=10, units="in", res=200)
  print(sda_hr_plot)
  dev.off()

  lst <- lapply(SDAdata, unlist)
  colnames(SDAdata_stepIntegral)<-colnames(SDAdata)
  SDAdata <- (data.frame(lapply(lst, `length<-`, max(lengths(lst)))))
  SDAdata<-rbind(SDAdata_stepIntegral, SDAdata)
  write.csv(file = SDAdata_name, SDAdata, row.names=FALSE)

}




