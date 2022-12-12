#' Title
#'
#' @param data.MMR The name of the MMR data file (“…analyzed.csv”; a character string); an output file from the MMR function.
#' @param data.SMR The name of the SMR data file (“…analyzed.csv”; a character string); an output file from the SMR function.
#' @param AnimalID Indicates individual ID; must be a vector of 4 characters. When missing, enter "NA"
#' @param BW.animal Indicates individual mass; must be a vector of 4 characters. When missing, enter "0"
#' @param resp.V Indicated the volume (L) of respirometer chambers; must be a vector of 4 numbers (e.g., c(1, 1, 1, 1), for four 1-L respirometers)
#' @param r2_threshold_smr R2 threshold for SMR, measurements below the threshold are excluded
#' @param r2_threshold_mmr R2 threshold for MMR, measurements below the threshold are excluded
#' @param min_length_mmr The duration of MMR steepest slope measurement; 180, 120, 90, 60 seconds (s)
#' @param scaling_exponent_mmr Body mass scaling exponent to correct MMR values for body size. MR=aBM^b (MR = metabolic rate, BW = body mass, a = scaling coefficient [the intercept], and b = scaling exponent [the power term])
#' @param scaling_exponent_smr Body mass scaling exponent to correct SMR values for body size. MR=aBM^b (MR = metabolic rate, BW = body mass, a = scaling coefficient [the intercept], and b = scaling exponent [the power term])
#' @param common_mass Metabolic performances are often calculated per unit mass. Use this argument to define what standardized mass should be. (default is MO2 mgO2kg-1^ min-1^, a common mass of 1 kg). Units = kg
#' @param mo2_val_for_calc Units of metabolic rates represented on export figures. mo2_common_mass_kg (standardized MO2 to 1 kg fish using scaling exponents provided); mg O2 min-1 kg-1), mo2_individual_kg (MO2 of the whole individual, mgO2 min-1)
#' @param plot_smr_quantile Indicating what percentile lower MO2 values to plot; must be 10, 15, or 20
#' @param date_format The date format used in the original FireSting data files. Must specify one of the following: "m/d/y", "d/m/y", "y-m-d"
#' @param N_Ch The number of channels of the FireSting. It must be either 4 or 8.
#' @param drop_ch Indicates which channel is dropped or entirely excluded from the analysis. Must be a numerical vector, e.g., c(1,3)
#' @param MLND Logical argument. If TRUE, SMR is estimated also using Mean Lowest Normal Distribution analysis. More details in Chabot et al 2016.
#' @param epoc_threshold Indicates the threshold relevant to an individual’s SMR to calculate the end time of recovery (the time at which metabolic rate has returned to epoc_threshold). The default is the SMR level (1; At 100 percent SMR). To use 120 percent SMR as a recovery threshold, enter epoc_threshold = 1.2.
#' @param recovMMR_threshold Short-term recovery threshold relevant to individual’s MMR. Indicates the (percent) level of MMR to calculate recovery indices. The default indicates 50 percent MMR, e.g., the time it takes to recover to 50 percent MMR (default). Enter 0.3 for 30 percent MMR, 0.8 for 80 percent MMR)
#' @param end_EPOC_Ch Manually assigned end time (min) of full recovery (EPOC). The time can be assigned for each channel independently; use a numerical vector of 4 variables, one for each channel (e.g., c(120, 120, NA, 180), for 2 h, 2h, SMR level, and 3 h EPOC end times, respectively)
#' @param background_prior The name of the analyzed background “…analyzed.csv” file, an output file from the SMR function (respiration in an empty respirometer measured before the respirometry trial).
#' @param background_post The same as background_prior, only post respirometry trial
#' @param background_slope Manually assigned background slope used to correct metabolic rate in all individuals. Provide numeric value in mgO2 L-1 h-1v
#' @param background.V Manually assigned respirometer volumes (L). A vector with 4 numeric variables, one for each channel.
#' @param background_linear_gr Logical argument. If TRUE, it is be assumed that bacterial growth (thus respiration rates) changed linearly across the duration of the respirometry trial. When TRUE, metabolic rate values across the given trial are corrected using the estimated background values from the linear growth curve. Both background_prior and background_post must be provided to apply this.
#' @param match_background_Ch Logical. If TRUE, the background respiration is estimated and applied channel-specific. The background_prior and background_post are used to estimate background respiration specific to each channel (or individual), which then is used to correct each individual’s MO2 independently. By default, the mean background respiration rate from all channels is calculated and applied to correct all individual’s MO2
#' @param mmr_background Specifies what background value should be used to correct MMR value. Options: i) "back_prior" takes background respiration rate value estimated from the background_prior file only, ii)"SAME_slope" takes background respiration rate value that is applied to the entire trial (MMR and SMR, indistinguishable), and iii) a defined respiration rate value mgO2 L-1 h-1.
#' @param local_path Logical. If TRUE (default) all returned files will be saved in the local working directory.
#' @param verbose.MLND From MLND: A logical controlling if a text progress bar from MLND is displayed during the fitting procedure. (see 'verbose' in mclust package functions).
#'
#' @importFrom stats lm coef var integrate predict quantile sd smooth.spline
#' @import graphics
#' @import grDevices
#' @importFrom gridExtra grid.arrange
#' @import scales
#' @import ggplot2
#' @import utils
#' @importFrom dplyr filter top_n arrange
#' @importFrom tidyr spread
#' @importFrom mclust Mclust densityMclust MclustBootstrap
#'
#' @return The output from \code{\link{print}}
#' @export
#'
MMR_SMR_AS_EPOC<-function(data.MMR = NULL,
                          data.SMR = NULL,
                          AnimalID,
                          BW.animal,
                          resp.V,
                          r2_threshold_smr,
                          r2_threshold_mmr,
                          min_length_mmr,
                          scaling_exponent_mmr = 1,
                          scaling_exponent_smr = 1,
                          common_mass = 1,
                          mo2_val_for_calc = "mo2_1kg",
                          plot_smr_quantile=15,
                          date_format = c("%Y-%m-%d %H:%M:%S", "GMT"),
                          N_Ch = 4,
                          drop_ch = NULL,
                          MLND = TRUE,
                          verbose.MLND = FALSE,
                          epoc_threshold = 1,
                          recovMMR_threshold = 0.5,
                          end_EPOC_Ch = c(NA, NA, NA, NA),
                          background_prior = NULL,
                          background_post = NULL,
                          background_slope = NULL,
                          background.V = NULL,
                          background_linear_gr = FALSE,
                          match_background_Ch = FALSE,
                          mmr_background = "SAME_slope",
                          local_path = TRUE){


  #  binding global variables locally to the function.
  DateTime_start<-m<-cycle_mmr<-cycle_type<-back_m_prior1<-back_m_prior2<-back_m_prior3<-back_m_prior4<-back_regression1<-back_regression2<-back_regression3<-back_regression4<-xpos<-ypos<-hjustvar<-vjustvar<-annotateText<-min_start<-quantiles<-mo2_perc<-smr_val<-smr_method<-NULL

  if(!length(as.vector(date_format))==2){
    stop_function<-TRUE
    if(stop_function) {
      stop("Argument 'date_format' is not properly stated. \n It must be a vector of two stating: i) the date time format and ii) timezone. \n Default is: c(\"%Y-%m-%d %H:%M:%S\", \"GMT\"). Argument is passed to strptime() ")
    }
  }

  ## recovery smoothing, only relevant withing this function
  ## # start of EPOC.spar plot function
  EPOC.spar <- function(spar,
                        d,
                        EPOCdata,
                        mmr.val,
                        epoc_threshold ,
                        recovMMR_threshold ,
                        newdata.smr,
                        MLND,
                        end_EPOC,
                        scaling_exponent_mmr,
                        scaling_exponent_smr,
                        common_mass){

  	#1 smooth and then predict the curve f
    fit<-smooth.spline(d$time_mo2,d$mo2, spar=spar)
  	f = function(x) {predict(fit, x)$y}

  	end<-round(d$time_mo2[nrow(d)],1)
  	newx<-seq(0,end, by=1)
  	newy<-predict(fit, newx, deriv=0)
  	lapply(newy, as.numeric)
  	newVal<-data.frame(Reduce(cbind,newy)) # creating dummy dataset with predicted values

  	smr.row<-newdata.smr[which(as.character(newdata.smr$Ch)==as.character(d$Ch[1])),]
  	ID<-smr.row["ID"]

  	# smr_mean10minVal SMR_low10quant SMR_low15quant SMR_low20quant  smr_mlnd
  	b1.1<-as.numeric(round(smr.row["smr_mean10minVal"],2))
  	b1.2<-as.numeric(round(smr.row["SMR_low10quant"],2))
  	b1.3<-as.numeric(round(smr.row["SMR_low15quant"],2))
  	b1.4<-as.numeric(round(smr.row["SMR_low20quant"],2))

  	if(MLND == TRUE){ # need this argument because MLND gives another SMR thresholdm that with MLND = FALSE is not available.
    	b1.5<-as.numeric(round(smr.row["smr_mlnd"],2))
  	}else{
  	  b1.5<-0
  	}

  	smr_type_list<-c("smr_mean10minVal", "SMR_low10quant", "SMR_low15quant", "SMR_low20quant", "smr_mlnd")
  	b_list<-c(b1.1,b1.2,b1.3,b1.4,b1.5)

  	for (i in 1:5){

  		if(i == 1){
  			EPOCdata.temp<-matrix(ncol=26, nrow=0)
  			colnames(EPOCdata.temp) <- c("ID", "smr_type", "smr","spar", "EPOC_full", "end_EPOC_min", "SMR_intergral_full", "SMR_threshold", "EPOC_1hr", "MO2_1hr", "EPOC_2hr", "MO2_2hr", "EPOC_3hr", "MO2_3hr", "EPOC_4hr", "MO2_4hr",  "EPOC_5hr", "MO2_5hr", "end_EPOC.mmr", "EPOC_mmr", "MO2_mmr", "MMR", "MMR_percent",  "scaling_exponent_mmr",  "scaling_exponent_smr", "common_mass")
  		}


  	  # adjustable epoc threshold (default is just SMR values as calculated using different methods)
  		if (epoc_threshold == 1){
  		  	b <- b_list[i]

  		}else{

  		    b <- b_list[i] * epoc_threshold
  		}

  	  b.mmr <- recovMMR_threshold * mmr.val # recv thresholds gives a proportion (or percent) and the mmr is the MMR is the MO2 value

  		smr_type <- smr_type_list[i]
  		#2-1- establish SMR threshold m1 and b1 for SMR (b1= y intercept, m1=slope)
  		m0 <- 0 # slope

  		# smr function
  		f.smr <- function(x)(m0*x)+b

  		# % mmr recovery function; default 50%
  		f.mmr <- function(x)(m0*x)+b.mmr # smr function

  		# if manually provided EPOC threshold does not exist then find when EPOC ends
  		if (is.na(end_EPOC)){
    		# Calculate the area under the SMR for all types of SMR calculated in the MMR_SMR_analyze function
    		# 2-2 the EPOC cutoff in the smoothed function / this is used also for the SMR block
    		end_EPOC <- newVal$init[(which(round(newVal$V2, 3)<=b))[1]]

    		if(is.na(end_EPOC)){
    			end_EPOC<-newVal$init[nrow(newVal)]
    		}
  		}

  		# % mmr end EPOC
  		end_EPOC.mmr <- newVal$init[(which(round(newVal$V2, 3)<=b.mmr))[1]]
  		if(is.na(end_EPOC.mmr)){
  			end_EPOC.mmr<-newVal$init[nrow(newVal)]
  			# message(paste("The animal does not reach % MMR recovery threshold: ", recovMMR_threshold, sep=""))
  		}

  		#2-3 the Breakpoit cuttoff
  		#bp1<-round(d$bp[1],0) # smr - or the y intercept ### This is the one Emily settled on for data analyses

  		f.int_mmr = function(x) {integrate(f, lower=0, upper=end_EPOC.mmr)$value}
  		f.int.vec_mmr = Vectorize(f.int_mmr, vectorize.args='x')
  		# f.int.vec_mmr = Vectorize(f.int_mmr, vectorize.args='x')
  		full_mmr<-f.int.vec_mmr(end_EPOC.mmr)
  		# SMR block
  		recovMMR<-integrate(f.mmr, lower=0, upper=end_EPOC.mmr)$value
  		EPOC_mmr<-round(full_mmr-recovMMR,3)
  		MO2_mmr<-round(newVal$V2[which(newVal$init==end_EPOC.mmr)],3)

  		#3 integrate recovery curve with to the provided end EPOC time
  		f.int = function(x) {integrate(f, lower=0, upper=end_EPOC)$value}
  		f.int.vec = Vectorize(f.int, vectorize.args='x')
  		f.int.vec = Vectorize(f.int, vectorize.args='x')
  		full<-f.int.vec(end_EPOC)
  		# SMR block
  		SMR<-integrate(f.smr, lower=0, upper=end_EPOC)$value
  		EPOC_full<-round(full-SMR,3)
  		MO2_full<-round(newVal$V2[which(newVal$init==end_EPOC)],3)

  		#3.1 integrate revovery curve for the first hr
  		f.int_1hr = function(x) {integrate(f, lower=0, upper=60)$value}
  		f.int.vec_1hr = Vectorize(f.int_1hr, vectorize.args='x')
  		full_1hr<-f.int.vec_1hr(60)
  		# SMR block
  		SMR_1hr<-integrate(f.smr, lower=0, upper=60)$value
  		EPOC_1hr<-round(full_1hr-SMR_1hr,3)
  		MO2_1hr<-round(newVal$V2[60],3)

  		#3.2 integrate revovery curve for the first 2hrs
  		f.int_2hr = function(x) {integrate(f, lower=0, upper=120)$value}
  		f.int.vec_2hr = Vectorize(f.int_2hr, vectorize.args='x')
  		full_2hr<-f.int.vec_2hr(120)
  		# SMR block
  		SMR_2hr<-integrate(f.smr, lower=0, upper=120)$value
  		EPOC_2hr<-round(full_2hr-SMR_2hr,3)
  		MO2_2hr<-round(newVal$V2[120],3)

  		#3.3 integrate revovery curve for the first 3hrs
  		f.int_3hr = function(x) {integrate(f, lower=0, upper=180)$value}
  		f.int.vec_3hr = Vectorize(f.int_3hr, vectorize.args='x')
  		full_3hr<-f.int.vec_3hr(180)
  		# SMR block
  		SMR_3hr<-integrate(f.smr, lower=0, upper=180)$value
  		EPOC_3hr<-round(full_3hr-SMR_3hr,3)
  		MO2_3hr<-round(newVal$V2[180],3)

  		#3.4 integrate revovery curve for the first 4hrs
  		f.int_4hr = function(x) {integrate(f, lower=0, upper=240)$value}
  		f.int.vec_4hr = Vectorize(f.int_4hr, vectorize.args='x')
  		full_4hr<-f.int.vec_4hr(240)
  		# SMR block
  		SMR_4hr<-integrate(f.smr, lower=0, upper=240)$value
  		EPOC_4hr<-round(full_4hr-SMR_4hr,3)
  		MO2_4hr<-round(newVal$V2[240],3)

  		#3.5 integrate revovery curve for the first 5hrs
  		f.int_5hr = function(x) {integrate(f, lower=0, upper=300)$value}
  		f.int.vec_5hr = Vectorize(f.int_5hr, vectorize.args='x')
  		full_5hr<-f.int.vec_5hr(300)
  		# SMR block
  		SMR_5hr<-integrate(f.smr, lower=0, upper=300)$value
  		EPOC_5hr<-round(full_5hr-SMR_5hr,3)
  		MO2_5hr<-round(newVal$V2[300],3)

  		# 3.bp1 integrate revovery curve for the first breakpoint
  		#f.int_bp1 = function(x) {integrate(f, lower=0, upper=bp1)$value}
  		#f.int.vec_bp1 = Vectorize(f.int_bp1, vectorize.args='x')
  		#full_bp1<-f.int.vec_bp1(bp1_val)
  		# SMR block
  		#SMR_bp1<-integrate(f.smr, lower=0, upper=bp1)$value
  		#EPOC_bp1<-round(full_bp1-SMR_bp1,2)
  		#MO2_bp1<-round(newVal$V2[bp1],2)
    	if(end_EPOC.mmr == 0){
    		  MO2_mmr=0
    	}

  		values<-as.data.frame(t(c(as.character(d$ID[1]),smr_type, b,  spar, EPOC_full, end_EPOC, SMR, epoc_threshold, EPOC_1hr, MO2_1hr, EPOC_2hr, MO2_2hr, EPOC_3hr, MO2_3hr, EPOC_4hr, MO2_4hr, EPOC_5hr, MO2_5hr, end_EPOC.mmr, EPOC_mmr, MO2_mmr, mmr.val, (MO2_mmr/mmr.val*100), scaling_exponent_mmr, scaling_exponent_smr, common_mass)))

  		colnames(values)<-c("ID", "smr_type", "smr","spar", "EPOC_full", "end_EPOC_min", "SMR_intergral_full", "SMR_threshold",
  		"EPOC_1hr", "MO2_1hr", "EPOC_2hr", "MO2_2hr", "EPOC_3hr", "MO2_3hr", "EPOC_4hr", "MO2_4hr",  "EPOC_5hr", "MO2_5hr", "end_EPOC.mmr", "EPOC_mmr", "MO2_mmr", "MMR", "MMR_percent", "scaling_exponent_mmr",  "scaling_exponent_smr", "common_mass")


  		EPOCdata.temp<-rbind(EPOCdata.temp,values)

  		col_smr<-c("black",  "darkmagenta", "darkred", "darkorange","darkgreen")
  		scale<-c(0.02,0.07,0.12,0.17,0.22)

  		if(MLND == TRUE){ # need this argument because MLND gives another SMR thresholdm that with MLND = FALSE is not available.
  			if(i==1){
  				plot(x=range(newVal$init), xlim=c(0,max(d$time_mo2, na.rm=TRUE)+50), ylim=c(0,max(d$mo2, na.rm=TRUE)),type='n', ylab="MO2", xlab="time (min)", main=paste("Smoothness=",spar, sep=""))
  				points(d$time_mo2,d$mo2)
  				lines(newy[[1]], newy[[2]], col="blue")
  				abline(v=60, col="grey", lty=2)
  				abline(v=120, col="grey", lty=2)
  				abline(v=180, col="grey", lty=2)
  				abline(v=240, col="grey", lty=2)
  				abline(v=300, col="grey", lty=2)
  				abline(h=b, col=col_smr[i],lty=1, lwd=1)
  				abline(v=end_EPOC, col=col_smr[i], lty=1)
  				abline(v=end_EPOC.mmr, col="magenta", lty=2,lwd=2)
  				text(x=(max(d$time_mo2)+50)-(0.01*(max(d$time_mo2)+50)),y=(max(d$mo2))-(scale[i]*(max(d$mo2))), label=paste("EPOC=",EPOC_full,"/ ",smr_type, sep=""), cex=0.8, col=col_smr[i], pos=2)

  			}else{
  				abline(h=b, col=col_smr[i],lty=1, lwd=1)
  				abline(v=end_EPOC, col=col_smr[i], lty=1)
  				abline(v=end_EPOC.mmr, col="magenta", lty=2,lwd=2)
  				text(x=(max(d$time_mo2)+50)-(0.01*(max(d$time_mo2)+50)),y=(max(d$mo2))-(scale[i]*(max(d$mo2))), label=paste("EPOC=",EPOC_full,"/ ",smr_type, sep=""), cex=0.8, col=col_smr[i], pos=2)

  			}

  		}else{
  			if(i==1){
  				plot(x=range(newVal$init), xlim=c(0,max(d$time_mo2, na.rm=TRUE)+50), ylim=c(0,max(d$mo2, na.rm=TRUE)),type='n', ylab="MO2", xlab="time (min)", main=paste("Smoothness=",spar, sep=""))
  				points(d$time_mo2,d$mo2)
  				lines(newy[[1]], newy[[2]], col="blue")
  				abline(v=60, col="grey", lty=2)
  				abline(v=120, col="grey", lty=2)
  				abline(v=180, col="grey", lty=2)
  				abline(v=240, col="grey", lty=2)
  				abline(v=300, col="grey", lty=2)
  				abline(h=b, col=col_smr[i],lty=1, lwd=1)
  				abline(v=end_EPOC, col=col_smr[i], lty=1)
  				abline(v=end_EPOC.mmr, col="magenta", lty=2,lwd=2)
  				text(x=(max(d$time_mo2)+50)-(0.01*(max(d$time_mo2)+50)),y=(max(d$mo2))-(scale[i]*(max(d$mo2))), label=paste("EPOC=",EPOC_full,"/ ",smr_type, sep=""), cex=0.8, col=col_smr[i], pos=2)


  			}
  		 if(!i==1 & !i==5){
  				abline(h=b, col=col_smr[i],lty=1, lwd=1)
  				abline(v=end_EPOC, col=col_smr[i], lty=1)
  				abline(v=end_EPOC.mmr, col="magenta", lty=2,lwd=2)
  				text(x=(max(d$time_mo2)+50)-(0.01*(max(d$time_mo2)+50)),y=(max(d$mo2))-(scale[i]*(max(d$mo2))), label=paste("EPOC=",EPOC_full,"/ ",smr_type, sep=""), cex=0.8, col=col_smr[i], pos=2)

  			}


  		}# end of "if(MLND == TRUE)"

  		if (i == 5){

  		  EPOCdata<-rbind(EPOCdata,EPOCdata.temp)
      }

  	}

  	return(EPOCdata)

}# end of EPOC.spar plot function



  graphics.off()

   filename.SMR<-paste(gsub('.{4}$', '',data.SMR), "_SMR", sep="")
	 filename.MMR<-paste(gsub('.{4}$', '',data.MMR), "_MMR_", sep="")


  # **********************************************
  # START-- >>> background
  if((is.null(background.V) & !is.null(background_slope)) |(!is.null(background.V) & is.null(background_slope))) {
    stop_function <- TRUE
    if(stop_function) stop("Must provide background unique slope and volume together, or provide a background file with same volumes as animal respirometers")
  }

	if(c(is.null(background_prior) & is.null(background_post)) & match_background_Ch){
    stop_function <- TRUE
	  if(stop_function) stop("If match_background_Ch = TRUE, must provide at least one datafiles, either background measured before or after the respirometry trial")
  }

  # background file wrangling
  # 1. find what channels recorded background
  if (!is.null(background_post) | !is.null(background_prior) ) {
    # 2 calculate mean for each background channel
    if((is.null(background.V) & !is.null(background_slope) ) |(!is.null(background.V) & is.null(background_slope))) {
      stop_function <- TRUE
      if(stop_function) stop("If using manually input background slopes, must provide both, a unique slope and volume of the repirometer")
    }

    if(!is.null(background_prior)){

      if(file.exists(background_prior) | file.exists(paste("./MMR_SMR_AS_EPOC/csv_input_files/", background_prior, sep=""))){ # after running through RMRrepeat - this will be saved in csv input files
      	  if(file.exists(paste("./MMR_SMR_AS_EPOC/csv_input_files/", background_prior, sep=""))){
            back_prior<-read.csv(paste("./MMR_SMR_AS_EPOC/csv_input_files/", background_prior, sep=""))
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

      # back_prior<-read.csv(background_prior)
      back_ch<-length(unique(back_prior$Ch))
    }


    if(!is.null(background_post)){
      if(file.exists(background_post) | file.exists(paste("./MMR_SMR_AS_EPOC/csv_input_files/", background_post, sep=""))){ # after running through RMRrepeat - this will be saved in csv input files
    	  if(file.exists(paste("./MMR_SMR_AS_EPOC/csv_input_files/", background_post, sep=""))){
          back_post<-read.csv(paste("./MMR_SMR_AS_EPOC/csv_input_files/", background_post, sep=""))
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

      back_ch<-length(unique(back_post$Ch))
    }

     # Jan 4 2020: make linear regression over time and correct background based on a predicted value
    # create
    if (background_linear_gr==TRUE){

      if(is.null(background_post) | is.null(background_post)){
       stop("Missing a file to model the growth of bacteria, must provide both: \"background_prior\" \"and background_post\" files")
      }

      if (!dir.exists("plots_background") | !dir.exists("./MMR_SMR_AS_EPOC/plots_background")){
        if(local_path & !dir.exists("plots_background")){
          dir.create(file.path("./plots_background"), recursive = TRUE)
        }
        if(local_path == FALSE &!dir.exists("./MMR_SMR_AS_EPOC/plots_background")){
          dir.create(file.path("./MMR_SMR_AS_EPOC/plots_background"), recursive = TRUE)
        }
      }

      back_prior$DateTime_start<- strptime(back_prior$DateTime_start, format = date_format[1], tz = date_format[2])
  	  back_post$DateTime_start<- strptime(back_post$DateTime_start, format = date_format[1], tz = date_format[2])

  	  back_all<- rbind(back_prior, back_post)
  	  back_all$DateTime_start<- as.POSIXct(back_all$DateTime_start)

  	  back_regression_plot<-ggplot(data=back_all, aes(x=DateTime_start, y=m, colour=Ch, group=Ch))+
        geom_point()+
        geom_smooth(method="lm", se=FALSE)+
        theme_bw()+
        ylab("Regression slope value- mgO2/min")+
  	    theme(axis.text.x = element_text(angle = 45))+
  	    facet_grid(Ch~.)

  	  if (local_path | !dir.exists("./MMR_SMR_AS_EPOC/plots_background")){
  	    plotname.backgr<-paste( filename.SMR,"_PLOT_BACKGROUND.png", sep="")
  	  }else{
  	   	plotname.backgr<-paste("./MMR_SMR_AS_EPOC/plots_background/", filename.SMR,"_PLOT_BACKGROUND.png", sep="")
  	  }

  	 if (match_background_Ch==TRUE){

  	    message("Background: assumed a linear change of bacterial respiration rates over the time of trial. Regressions are specific to the channel (respirometer)")

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
  	    message("Background: assumed a linear change of bacterial respiration rates over the time of trial. One regression for all channels")
      }


    } # end for getting linear regressions for the background

    if (background_linear_gr==FALSE){
      if (!is.null(background_prior)){
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
      # 3. estimate one background slope mean to be used in MR corrections
      if (!is.null(background_post)){

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

        if (!is.null(background_prior)){
          ch_available<-as.numeric(substr(as.character(unique(back_prior$Ch)),start=3, stop=3))
        }
        if (!is.null(background_post)){
          ch_available<-as.numeric(substr(as.character(unique(back_post$Ch)),start=3, stop=3))
        }


        for (i in 1:back_ch){

          # prior only condition
         if(exists(paste("back_m_prior", ch_available[i], sep="")) & !exists(paste("back_m_post", ch_available[i], sep=""))){

           if(back_ch_prior_names[[i]]==paste("back_m_prior", ch_available[i], sep="")){
              # message("matching background Channels")

              back_m_name2<-paste("back_m", ch_available[i], sep="")
              assign(back_m_name2, back_ch_prior[[i]])
              #next # continue the loop to the nect channel

           }

         }

          # post only condition
          if(!exists(paste("back_m_prior", ch_available[i], sep="")) & exists(paste("back_m_post",ch_available[i], sep=""))){

            if(back_ch_post_names[[i]]==paste("back_m_post", ch_available[i], sep="")){
              # message("matching background Channels")

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

        if(is.null(background_post) | !is.null(background_prior)){

          prior_mean<-sum(as.numeric(back_ch_prior)) / (length (unique(back_prior$Ch)))
          back_m<-prior_mean

        }

        if(!is.null(background_post) | is.null(background_prior)){

          post_mean<-sum(as.numeric(back_ch_post)) / (length (unique(back_post$Ch)))
          back_m<-post_mean

        }

        if(!is.null(background_post) & !is.null(background_prior)){

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
  # **********************************************



	newdata.smr<-as.data.frame(matrix(ncol=33, nrow=0))
	names(newdata.smr)<-c("filename", "ID", "Ch", "BW","t_min","t_max", "t_mean", "N_mo2", #8
	"smr_mean10minVal","smr_SD10minVal", "smr_CV10minVal", "SMR_low10quant","SMR_low15quant","SMR_low20quant", #6
	"smr_mlnd", "smr_CVmlnd", "smr_Nmlnd", #3
	"mmr", "mmr_overall", #2
	"AS_smr_mean10minVal", "AS_SMR_low10quant", "AS_SMR_low15quant", "AS_SMR_low20quant", "AS_smr_mlnd", #5
	"AS_smr_mean10minVal_overall", "AS_SMR_low10quant_overall", "AS_SMR_low15quant_overall", "AS_SMR_low20quant_overall", "AS_smr_mlnd_overall", #5
	"mmr_length_cycle", "scaling_exponent_mmr", "scaling_exponent_smr", "common_mass")#1

	# cols = c(4:17, 19:32)
	# newdata.smr[,cols] %<>% lapply(function(x) as.numeric(as.character(x)))
	# cols2 = c(1:3, 18)
	# newdata.smr[,cols2] %<>% lapply(function(x) as.character(x))

	cols = c(4:33)
	cols2 = c(1:3)
  newdata.smr[,cols] <- lapply(newdata.smr[,cols], as.character)
	newdata.smr[,cols] <- lapply(newdata.smr[,cols], as.numeric)
	newdata.smr[,cols2] <- lapply(newdata.smr[,cols2], as.character)

 	# **********************************************
  # START -- >>> MMR file available
  if (!is.null(data.MMR)){

    if(file.exists(data.MMR) | file.exists(paste("./MMR_SMR_AS_EPOC/csv_input_files/", data.MMR, sep=""))){ # after running through RMRrepeat - this will be saved in csv input files
    	  if(file.exists(paste("./MMR_SMR_AS_EPOC/csv_input_files/", data.MMR, sep=""))){
          d_MMR<-read.csv(paste("./MMR_SMR_AS_EPOC/csv_input_files/", data.MMR, sep=""))
        }
        if(file.exists(data.MMR)){
          d_MMR<-read.csv(data.MMR)
        }
    	}else{
        stop_function<-TRUE
        if(stop_function){
          stop("Cannot locate the indicated data.MMR data file.")
        }
    }
    # d_MMR<-read.csv(data.MMR)

    # feb 23 note:
    # the previous version of the code before MMR became a flexible version (beginning of feb 2020), has an extra column: the delay is still there. If that is still there, but that file is desired to be used in the "new" MMR_SMR_AS_EPOC, thne get rid of that column and write a warning message:
    if(colnames(d_MMR)[4] == "delay_min"){
      d_MMR<-d_MMR[, -4]
      message("Note: an older version of MMR analyzed file was used: dropping unused column \"delay_min\"")
    }


    ### getting the right data frame of MMR values / choosing the right "sliding"
  	if (!c(min_length_mmr == 1 | min_length_mmr == 60 | min_length_mmr == 90 | min_length_mmr == 120 | min_length_mmr == 180)){
      message(paste("MMR: not identified duration of MMR measurement:", min_length_mmr, "seconds, re-assigning it to min_length_mmr = 1 (full measurement cycle)"))
  	  min_length_mmr<-1
  	 }

  	times_list<-c(60,90,120,180,1)
  	ntime<-which(times_list==min_length_mmr)

    d_MMR$DateTime_start<- strptime(d_MMR$DateTime_start, format = date_format[1], tz = date_format[2])

  	# drop any unwanted channels
    	if(!is.null(drop_ch[1])){
    	    n_ch_drop<-length(drop_ch)
    	     for(i in 1:n_ch_drop){
    	      d_MMR<-d_MMR[(substr(as.character(d_MMR$Ch), start=3, stop=3)!=drop_ch[i]),]
    	      }
    	}

  	d_MMR$cycle_type<-factor(d_MMR$cycle_type)
  	d_MMR$Ch<-factor(d_MMR$Ch)


  	d_MMR$cycle_mmr[d_MMR$cycle_mmr=="1" | d_MMR$cycle_mmr=="2"| d_MMR$cycle_mmr=="3"| d_MMR$cycle_mmr=="4"]<-1 # 1 is an arbritrary number indicating that this is a "full length MMR file slope"
  	d_MMR$cycle_mmr[(d_MMR$cycle_mmr=="1" | d_MMR$cycle_mmr=="2"| d_MMR$cycle_mmr=="3"| d_MMR$cycle_mmr=="4") & grepl("cycle", as.character(d_MMR$cycle_type))]<-10 # 10 is an arbitrary number that is not going to show up anywhere else/ needed to indicate teh cycle "full length" slope

  	#### SELECT the type of MMR SLOPES as wanted
  	# select either the absolute steepest slopes or the mean lowest slope
  	# See Zhang et al PEAk slopes JEB paper -- incorporate
    	# if (mmr_type == "min"){
    	#   d_MMR<-d_MMR[d_MMR$cycle_type!="Mean_1minMean", ]
    	#   levels(droplevels(d_MMR$cycle_type))
    	# }
  	  # }else{ # this will be the autocorrelation slope
  	  # d_MMR<-d_MMR[d_MMR$cycle_type!="Mean_1minSlopes", ]
  	  # levels(droplevels(d_MMR$cycle_type))
  	  # }
  	###

  	d<-d_MMR %>%
  		group_by(Ch) %>%
  		filter(cycle_mmr== min_length_mmr| cycle_mmr== 1 | cycle_type=="cycle2" | cycle_type=="cycle3" | cycle_type=="cycle4")
  	d<-as.data.frame(d)

  	d_MMR_list<-split(d_MMR, d_MMR$Ch)

  	d_list<-split(d, d$Ch)

    d_temp1<-as.data.frame(matrix(nrow=0, ncol=13))
    colnames(d_temp1)<-c("cycle_type", "cycle_start","cycle_end",  "cycle_mmr", "r2" ,"m", "b" , "t_min", "t_max", "t_mean", "Ch", "DateTime_start")

  	for (i in 1:length(d_list)){
  	  d_ch<-as.data.frame(d_list[i])
  	  colnames(d_ch)<-c("cycle_type", "cycle_start","cycle_end", "cycle_mmr", "r2" ,"m", "b" , "t_min", "t_max", "t_mean", "Ch", "DateTime_start")

  	  # d_temp<-d_ch
  	  # d_temp<-d_ch[grepl("MMR", as.character(d_ch$cycle_type)),]
  	  if((nrow(as.data.frame(d_ch[grepl("MMR", as.character(d_ch$cycle_type)),])))>1){

  	    d_ch<-d_ch[c(which((as.numeric(d_ch$cycle_mmr)>1 & grepl("MMR", as.character(d_ch$cycle_type))) | grepl("cycle", as.character(d_ch$cycle_type)))),]
  	  }
  	  d_temp1<-rbind(d_temp1, d_ch)

  	}

    d<-d_temp1

  	if(any(d$r2<r2_threshold_mmr)){
  		d_new<-d_MMR[0,]

  		for (i in 1:length(d_MMR_list)){
  			d_ch<-as.data.frame(d_MMR_list[i])
  			colnames(d_ch)<-c("cycle_type", "cycle_start","cycle_end", "cycle_mmr", "r2" ,"m", "b" , "t_min", "t_max", "t_mean", "Ch", "DateTime_start")

  			d_temp<-d_ch[grepl("MMR", as.character(d_ch$cycle_type)),]

  			if (nrow(d_temp) == 1 & c(d_temp$cycle_end[1] - d_temp$cycle_start[1] < 3) &
  			    c(nrow(d_temp) == 1)){ # no data with r2 at the threshold and at the desired minute length
            message(paste("MMR: ", d_temp$Ch[1], ": measurement was < 3 minutes, use full measurement cycle", sep =""))
  			}

  			if (all(d_temp$r2 < r2_threshold_mmr) | c(nrow(d_temp) == 1 & c(d_temp$cycle_end[1] - d_temp$cycle_start[1] < 3))){

  				if (all(d_temp$r2 < r2_threshold_mmr)){
  						new_min_length_mmr<-d_temp$cycle_mmr[d_temp$r2==max(d_temp$r2)][1]
  						d_ch<-d_ch[c(which(grepl("cycle", as.character(d_ch$cycle_type)) | d_ch$cycle_mmr== new_min_length_mmr)),]
  						d_new<-rbind(d_new, d_ch)

  					if(new_min_length_mmr == 1){
						  message(paste("MMR: ", d_ch$Ch[1],": chosen MMR r2 is too high, no data with: ", r2_threshold_mmr," | insted used the full measurement cycle with highest r2 = ", d_ch$r2[1], sep=""))
						}else{
						  message(paste("MMR: ", d_ch$Ch[1],": chosen MMR r2 is too high, no data with: ", r2_threshold_mmr," | instead used the", new_min_length_mmr ," s measurement cycle with highest r2 = ", d_ch$r2[1], sep=""))
						}
  						# message(paste(d_ch$Ch[1],": MMR measures all r2 below the set threshold: ", r2_threshold_mmr," / USE max r2 = ", d_ch$r2[1], sep=""))
  						# message(paste(d_ch$Ch[1],": MMR measure time extended from ", min_length_mmr ," to ", new_min_length_mmr , sep=""))
  				}

  			  if(!nrow(d_temp) == 1){
  			    for(j in ntime:(length(times_list))){
    					if(d_temp$r2[d_temp$cycle_mm==times_list[j]]>=r2_threshold_mmr){
    						new_min_length_mmr<-times_list[j]
    						d_ch<-d_ch[c(which(grepl("cycle", as.character(d_ch$cycle_type)) | d_ch$cycle_mmr == new_min_length_mmr)),]
    						d_new<-rbind(d_new, d_ch)

    							if (new_min_length_mmr==1){
    								message(paste("MMR: ", d_ch$Ch[1],": measurement time extended from ", min_length_mmr ,"s to full a measurement", sep=""))
    							}else{
    								message(paste("MMR: ", d_ch$Ch[1],": measurement time extended from ", min_length_mmr ,"s to ", new_min_length_mmr ,"s MMR measurement window", sep=""))
    							}

    						break
    					}
    				}
  			  }


  			}else{
  				d_ch<-d_ch[c(which(grepl("cycle", as.character(d_ch$cycle_type)) | d_ch$cycle_mmr==min_length_mmr)),]
  				d_new<-rbind(d_new, d_ch)
  			}
  		}

  		d_new$cycle_type[d_new$cycle_type=="MMR_slide"]<-"MMR"

  	}else{

  		d$cycle_type[d$cycle_type=="MMR_slide"]<-"MMR"
  		d_new<-d
  	}

  	d_MMR<-as.data.frame(d_new)
  	### MMR data set established

  	if (!is.null(data.MMR)){
  		d_MMR$bw<-NA
  		d_MMR$mo2<-NA
  		d_MMR$ID<-NA
      d_MMR$resp.V<-NA
      d_MMR$scaling_exponent<-scaling_exponent_mmr
      d_MMR$common_mass<-common_mass
  	}

  	for(i in 1:4){
			if( any(grepl(as.character(i),as.character(d_MMR$Ch)))){
				nameCh<-paste("Ch",i,sep="")
				bw.val<-BW.animal[i]
				ID<-AnimalID[i]
				resp.Vol<-resp.V[i]
				n.row<-which(d_MMR$Ch==nameCh)
				d_MMR$bw[n.row]<-bw.val
				d_MMR$ID[n.row]<-as.character(ID)
				d_MMR$resp.V[n.row]<-resp.Vol
			}
		}



  	# **********************************************
    # START -- >>> background corrections MMR
  	# no change with the linear background additions. MMR is the same since it is one measurement at the beginning.
  	# MMR end values MR in mgO2/min/kg


    # mmr_background = "SAME_slope" (default is to take the same value as the entire smr), must specify mmr_background = "back_prior" to have a prior slope only, or add a specific value to be substracted for mmr background

  	 d_MMR$DateTime_start<- strptime(d_MMR$DateTime_start, format="%Y-%m-%d %H:%M:%S")

     if(mmr_background != "SAME_slope" & background_linear_gr==FALSE){

        # message("MMR slope corrected using the same background slope value and methods as for all SMR slopes")

        if(mmr_background == "back_prior"){

          if(match_background_Ch != TRUE){ # if NOT channel specific

            if(!exists("prior_mean")) stop("Incomplete request: no prior-trail background respiration data to correct MMR measurement, must provide \"background_prior\"file")

              for (i in 1:nrow(d_MMR)){
          		  d_MMR$mo2[i]<-(( d_MMR$m[i]*(d_MMR$resp.V[i]-d_MMR$bw[i])) - (prior_mean * d_MMR$resp.V[i])) /(d_MMR$bw[i])#^scaling_exponent_mmr) # units mgO2 kg-1 min-1 -CORRECTED for back resp
              }

            message("MMR corrected for background: using a common background from a \"background_prior\" file only (same slope for all channels) | bacterial respiration slope: ", prior_mean)

            }else{ # this will be mmr slope different than smr, and is channel specific, i.e. match__background==TRUE
              # back_m_prior1
              # back_m_prior2
              # back_m_prior3
              # back_m_prior4
             message("MMR corrected for background: using Ch specific average background values from \"background_prior\"")
               for (i in 1:nrow(d_MMR)){

        		    if(substr(d_MMR$Ch[i], start=3, stop=3) == "1"){
        		    d_MMR$mo2[i]<-(( d_MMR$m[i]*(d_MMR$resp.V[i]-d_MMR$bw[i])) - (back_m_prior1 * d_MMR$resp.V[i])) /(d_MMR$bw[i])#^scaling_exponent_mmr) # units mgO2 kg-1 min-1 -CORRECTED for back resp
        		    # message("Channel: 1,  bacterial respiration slope: ", back_m1)
        		    }

        		    if(substr(d_MMR$Ch[i], start=3, stop=3) == "2"){
        		      d_MMR$mo2[i]<-(( d_MMR$m[i]*(d_MMR$resp.V[i]-d_MMR$bw[i])) - (back_m_prior2 * d_MMR$resp.V[i])) /(d_MMR$bw[i])#^scaling_exponent_mmr) # units mgO2 kg-1 min-1 -CORRECTED for back resp
        		      # message("Channel: 2,  bacterial respiration slope: ", back_m2)
        		    }

        		    if(substr(d_MMR$Ch[i], start=3, stop=3) == "3"){
        		      d_MMR$mo2[i]<-(( d_MMR$m[i]*(d_MMR$resp.V[i]-d_MMR$bw[i])) - (back_m_prior3 * d_MMR$resp.V[i])) /(d_MMR$bw[i])#^scaling_exponent_mmr) # units mgO2 kg-1 min-1 -CORRECTED for back resp
        		      # message("Channel: 3,  bacterial respiration slope: ", back_m3)
        		    }

        		    if(substr(d_MMR$Ch[i], start=3, stop=3) == "4"){
        		      d_MMR$mo2[i]<-(( d_MMR$m[i]*(d_MMR$resp.V[i]-d_MMR$bw[i])) - (back_m_prior4 * d_MMR$resp.V[i])) /(d_MMR$bw[i])#^scaling_exponent_mmr) # units mgO2 kg-1 min-1 -CORRECTED for back resp
        		      # message("Channel: 4,  bacterial respiration slope: ", back_m4)
        		    }

        		  }

          }

        }else{ # else from mmr = back prior

          if(!is.numeric(mmr_background)) stop("Incomplete request: must provide \"mmr_background\" slope value to manually MMR for background OR use mmr_background=\"SAME_slope\"")

           for (i in 1:nrow(d_MMR)){
        		  d_MMR$mo2[i]<-(( d_MMR$m[i]*(d_MMR$resp.V[i]-d_MMR$bw[i])) - (mmr_background * d_MMR$resp.V[i])) /(d_MMR$bw[i])#^scaling_exponent_mmr) # units mgO2 kg-1 min-1 -CORRECTED for back resp
           }
          message("MMR corrected for background: using specified slope value, same slope for all channels | bacterial respiration slope: ", mmr_background)
        }


		  }else{ # ALL OF THE CORRECTION BELOW IS WITH USING SAME SLOPE AS SMR

		    if(background_linear_gr==TRUE){
		        if(match_background_Ch==TRUE){
		          message("MMR corrected for background: using Ch specific background based on estimated linear regression of bacterial growth")
		            for (i in 1:nrow(d_MMR)){

		              background_slopes<-data.frame(matrix(ncol=1, nrow=nrow(d_MMR)))
		              colnames(background_slopes)<-c("back_m")

          		    if(substr(as.character(d_MMR$Ch[i]), start=3, stop=3) == "1"){
          		      background_slopes$back_m[i]<-predict(back_regression1, data.frame(DateTime_start = d_MMR$DateTime_start[i]))
                    back_m1<-predict(back_regression1, data.frame(DateTime_start = d_MMR$DateTime_start[i]))
            		    d_MMR$mo2[i]<-(( d_MMR$m[i]*(d_MMR$resp.V[i]-d_MMR$bw[i])) - (back_m1 * d_MMR$resp.V[i])) /(d_MMR$bw[i])#^scaling_exponent_mmr) # units mgO2 kg-1 min-1 -CORRECTED for back resp
          		      # message("Channel: 1,  bacterial respiration slope: ", back_m1)
          		    }

          		    if(substr(as.character(d_MMR$Ch[i]), start=3, stop=3) == "2"){
          		      background_slopes$back_m[i]<-predict(back_regression2, data.frame(DateTime_start = d_MMR$DateTime_start[i]))
          		      back_m2<-predict(back_regression2, data.frame(DateTime_start = d_MMR$DateTime_start[i]))
          		      d_MMR$mo2[i]<-(( d_MMR$m[i]*(d_MMR$resp.V[i]-d_MMR$bw[i])) - (back_m2 * d_MMR$resp.V[i])) /(d_MMR$bw[i])#^scaling_exponent_mmr) # units mgO2 kg-1 min-1 -CORRECTED for back resp
          		      # message("Channel: 2,  bacterial respiration slope: ", back_m2)
          		    }

          		    if(substr(as.character(d_MMR$Ch[i]), start=3, stop=3) == "3"){
          		      background_slopes$back_m[i]<-predict(back_regression3, data.frame(DateTime_start = d_MMR$DateTime_start[i]))
          		      back_m3<-predict(back_regression3, data.frame(DateTime_start = d_MMR$DateTime_start[i]))
          		      d_MMR$mo2[i]<-(( d_MMR$m[i]*(d_MMR$resp.V[i]-d_MMR$bw[i])) - (back_m3 * d_MMR$resp.V[i])) /(d_MMR$bw[i])#^scaling_exponent_mmr) # units mgO2 kg-1 min-1 -CORRECTED for back resp
          		      # message("Channel: 3,  bacterial respiration slope: ", back_m3)
          		    }

          		    if(substr(as.character(d_MMR$Ch[i]), start=3, stop=3) == "4"){
                    background_slopes$back_m[i]<-predict(back_regression4, data.frame(DateTime_start = d_MMR$DateTime_start[i]))
          		      back_m4<-predict(back_regression4, data.frame(DateTime_start = d_MMR$DateTime_start[i]))
          		      d_MMR$mo2[i]<-(( d_MMR$m[i]*(d_MMR$resp.V[i]-d_MMR$bw[i])) - (back_m4 * d_MMR$resp.V[i])) /(d_MMR$bw[i])#^scaling_exponent_mmr) # units mgO2 kg-1 min-1 -CORRECTED for back resp
          		      # message("Channel: 4,  bacterial respiration slope: ", back_m4)
          		    }
          		  }# end of the looop
		        }
            if(match_background_Ch==FALSE){
                background_slopes<-data.frame(DateTime_start = d_MMR$DateTime_start)
                background_slopes$back_m<-predict(back_regression, background_slopes)

            	 for (i in 1:nrow(d_MMR)){
            	  d_MMR$mo2[i]<-(( d_MMR$m[i]*(d_MMR$resp.V[i]-d_MMR$bw[i])) - (background_slopes$back_m[i] * d_MMR$resp.V[i])) /(d_MMR$bw[i])#^scaling_exponent_mmr) # units mgO2 kg-1 min-1 -CORRECTED for back resp
            	 }
		        }


		    }else{
            	# 1.1 if background files (either prior or post, or both) are provided and its one overall mean value (back_m)
        		if ((( !is.null(background_post) | !is.null(background_prior)) & match_background_Ch==FALSE) & is.null(background_slope) ){


        		  for (i in 1:nrow(d_MMR)){
        		    d_MMR$mo2[i]<-(( d_MMR$m[i]*(d_MMR$resp.V[i]-d_MMR$bw[i])) - (back_m * d_MMR$resp.V[i])) /(d_MMR$bw[i])#^scaling_exponent_mmr) # units mgO2 kg-1 min-1 -CORRECTED for back resp
        		  }
        		  # message("Correcting for average background from file - MMR corrected (same slope value as for SMR) ")
        		  message("MMR corrected for background: using a mean (prior and/or post) background measurements | mean bacterial respiration slope: ", back_m)

        		}

        		## if background slope and volume are specifically provided, then use those! this alos overrides the background prior and post argument.
        		# all channels with the same slope
        		if (!is.null(background_slope) ){

        		  for (i in 1:nrow(d_MMR)){
        		    d_MMR$mo2[i]<-(( d_MMR$m[i]*(d_MMR$resp.V[i]-d_MMR$bw[i])) - (background_slope * background.V)) /(d_MMR$bw[i])#^scaling_exponent_mmr) # units mgO2 kg-1 min-1 -CORRECTED for back resp
        		  }
        		  # message("Correcting for common background slope - only MMR corrected (same slope value as for SMR) ")
        		  message("MMR corrected for background: using a common manually provided background slope for all channels | bacterial respiration slope: ", background_slope)
        		}

        		# 1.2 if background files are provided and its channel specific - MEAN SLOPE OF PRIOR AND POST
        		if (( !is.null(background_post) | !is.null(background_prior)) & match_background_Ch==TRUE){
        		  message("MMR: correcting for Ch specific average background (using mean slope back prior and/or post), same slope value as for SMR")
               for (i in 1:nrow(d_MMR)){

        		    if(substr(d_MMR$Ch[i], start=3, stop=3) == "1"){
        		    d_MMR$mo2[i]<-(( d_MMR$m[i]*(d_MMR$resp.V[i]-d_MMR$bw[i])) - (back_m1 * d_MMR$resp.V[i])) /(d_MMR$bw[i])#^scaling_exponent_mmr) # units mgO2 kg-1 min-1 -CORRECTED for back resp
        		    # message("Channel: 1,  bacterial respiration slope: ", back_m1)

        		    }

        		    if(substr(d_MMR$Ch[i], start=3, stop=3) == "2"){
        		      d_MMR$mo2[i]<-(( d_MMR$m[i]*(d_MMR$resp.V[i]-d_MMR$bw[i])) - (back_m2 * d_MMR$resp.V[i])) /(d_MMR$bw[i])#^scaling_exponent_mmr) # units mgO2 kg-1 min-1 -CORRECTED for back resp
        		      # message("Channel: 2,  bacterial respiration slope: ", back_m2)

        		    }

        		    if(substr(d_MMR$Ch[i], start=3, stop=3) == "3"){
        		      d_MMR$mo2[i]<-(( d_MMR$m[i]*(d_MMR$resp.V[i]-d_MMR$bw[i])) - (back_m3 * d_MMR$resp.V[i])) /(d_MMR$bw[i])#^scaling_exponent_mmr) # units mgO2 kg-1 min-1 -CORRECTED for back resp
        		      # message("Channel: 3,  bacterial respiration slope: ", back_m3)

        		    }

        		    if(substr(d_MMR$Ch[i], start=3, stop=3) == "4"){
        		      d_MMR$mo2[i]<-(( d_MMR$m[i]*(d_MMR$resp.V[i]-d_MMR$bw[i])) - (back_m4 * d_MMR$resp.V[i])) /(d_MMR$bw[i])#^scaling_exponent_mmr) # units mgO2 kg-1 min-1 -CORRECTED for back resp
        		      # message("Channel: 4,  bacterial respiration slope: ", back_m4)

        		    }
        		  }
        		}


        		# 2. if background files are not provided
        		if ((is.null(background_post) & is.null(background_prior)) & is.null(background_slope)){
        		  message("MMR: no correction for background respiration")

          		for (i in 1:nrow(d_MMR)){
        		  	d_MMR$mo2[i]<-d_MMR$m[i]*(d_MMR$resp.V[i]-d_MMR$bw[i])/(d_MMR$bw[i])#^scaling_exponent_mmr) # units mgO2 kg-1 min-1
          		}
        		}
		    }
		  }

    # END -- >>> background corrections MMR
  	# **********************************************


  	# size correct each value for SMR calculations, lead with this into EPOC calculations
		d_MMR$mo2_individual_kg<-d_MMR$mo2 * d_MMR$bw
		# d_MMR$mo2_common_mass_kg<-exp( (scaling_exponent_mmr*(log(d_MMR$common_mass)-log(d_MMR$bw))) + log(d_MMR$mo2_individual_kg))
		d_MMR$mo2_common_mass_kg<-exp( (scaling_exponent_mmr*(log(d_MMR$common_mass)-log(d_MMR$bw))) + log(d_MMR$mo2_individual_kg))

    # ScalingSlope = (y2 - y1) / (x2 - x1)
    # y2 = log (mo2) <<< current mass dependent value
    # y1 = unknown <<< solve for; this will be log() values! use (exp()) to back-transform
    # x2 = log (mass kg) <<< current body mass in kg
    # x1 = log (common_mass) <<< common body mass for the correction

    # solve for
    # y1 = y2 - ( ScalingSlope * (x2 - x1) )
		# or solve for y2:
		# y2 = ScalingSlope * (x2 - x1) + y1
	  d_MMR$mo2_common_mass_kg <- exp( log(d_MMR$mo2_individual_kg) - (scaling_exponent_mmr * (log(d_MMR$bw) - log(d_MMR$common_mass)) ) )

	  # rename the mo2 to indicate that this is per kg fish
	  names(d_MMR)[names(d_MMR) == 'mo2'] <- 'mo2_1kg'

    d_MMR$mo2<-d_MMR[,mo2_val_for_calc]

  }
	# END -- >>> MMR file available
 	# **********************************************



	# **********************************************
	# START -- >>> If SMR data IS available
	# analyses MMR points for SMR value distributions and calculations too if available, no EPOC
	if(!is.null(data.SMR)){

    if(file.exists(data.SMR) | file.exists(paste("./MMR_SMR_AS_EPOC/csv_input_files/", data.SMR, sep=""))){ # after running through RMRrepeat - this will be saved in csv input files
    	  if(file.exists(paste("./MMR_SMR_AS_EPOC/csv_input_files/", data.SMR, sep=""))){
          d_SMR<-read.csv(paste("./MMR_SMR_AS_EPOC/csv_input_files/", data.SMR, sep=""))
        }
        if(file.exists(data.SMR)){
          d_SMR<-read.csv(data.SMR)
        }
    	}else{
        stop_function<-TRUE
        if(stop_function){
          stop("Cannot locate the indicated data.SMR data file.")
        }
    }


		# d_SMR<-read.csv(data.SMR)
			# drop any unwanted channels
	  if(!is.null(drop_ch[1])){
	    n_ch_drop<-length(drop_ch)
	     for(i in 1:n_ch_drop){
	      d_SMR<-d_SMR[(substr(as.character(d_SMR$Ch), start=3, stop=3)!=drop_ch[i]),]
	      }
	  }
    d_SMR$Ch<-factor(d_SMR$Ch)


  	if (local_path | !dir.exists("MMR_SMR_AS_EPOC")){
  		plotname.freq<-paste( filename.SMR,"_PLOT_SMR_analyses.png", sep="")
  		plotname.smr.meth<-paste( filename.SMR,"_PLOT_SMR_methodsALL.png", sep="")
  	}else{
  		plotname.freq<-paste("./MMR_SMR_AS_EPOC/plots_min_values_SMR/", filename.SMR,"_PLOT_SMR_analyses.png", sep="")
  		plotname.smr.meth<-paste("./MMR_SMR_AS_EPOC/plots_methods_sum_SMR/", filename.SMR,"_PLOT_SMR_methodsALL.png", sep="")
  	}

		if(!colnames(d_SMR)[11]=="type"){
			d_SMR$type="SMR"
			d_SMR<-d_SMR[,c("time_frame", "min_start", "r2", "b", "m", "t_min", "t_max", "t_mean" ,"Ch", "DateTime_start", "type", "n_min", "ID_code" )]
		}


		## choose what to keep and what not for the SMR
		# 2) keep only > 1 min sections for SMR calculations any type  selected above (SMR, pre-, post-shallow slopes) and exclude "SMR-cut", "SMR-cut1", "SMR-cut2"
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
		  stop_function<-TRUE
        if(stop_function){
          stop(paste("SMR/RMR: ! NO DATA: The r2 = ", (r2_threshold_smr), " is too high. The highest r2 for SMR/RMR measurement is ", min(d_SMR$r2), sep=""))
        }
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
		d_SMR$common_mass<-common_mass

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


    ## this is where lines indicated with NA in the individual animal argument are ditched
#     d_SMR<-d_SMR[d_SMR$ID!="NA",]
#     d_MMR<-d_MMR[d_MMR$ID!="NA",]
#
# 		d_MMR<-d_MMR[as.character(d_MMR$ID)!="blank",]
# 		d_SMR<-d_SMR[as.character(d_SMR$ID)!="blank",]
#
#     d_SMR$ID<-factor(d_SMR$ID)
#     d_MMR$ID<-factor(d_MMR$ID)




		# **********************************************
    # START -- >>> background corrections SMR
		# Jan 5 2020 change - account for background using linear regression
    # making predictions

# 		d_SMR$DateTime_start<-as.character(gsub("(", "", d_SMR$DateTime_start, fixed=TRUE))
#   	d_SMR$DateTime_start<-as.character(gsub(")", "", d_SMR$DateTime_start, fixed=TRUE))

	   d_SMR$DateTime_start<- strptime(d_SMR$DateTime_start, format = date_format[1], tz = date_format[2])

#   	if (date_format== "d/m/y"){
# 	     d_SMR$DateTime_start<- strptime(d_SMR$DateTime_start, format="%d/%m/%y %H:%M:%OS")
#   	}
#   	if (date_format== "y-m-d"){
#       # DateTime<- chron(dates=back_$date,times=back_$time,format=c('y-m-d','h:m:s'))
# 	     d_SMR$DateTime_start<- strptime(d_SMR$DateTime_start, format="%y-%m-%d %H:%M:%OS")
#   	}

		if(background_linear_gr==TRUE & match_background_Ch==FALSE){
        background_slopes<-data.frame(DateTime_start = d_SMR$DateTime_start)
        background_slopes$back_m<-predict(back_regression, background_slopes)
  	   for (i in 1:nrow(d_SMR)){
  	    d_SMR$mo2[i]<-(( d_SMR$m[i]*(d_SMR$resp.V[i]-d_SMR$bw[i])) - (background_slopes$back_m[i] * d_SMR$resp.V[i])) /(d_SMR$bw[i])#^scaling_exponent_smr) # units mgO2 kg-1 min-1 -CORRECTED for back resp
  	   }

      background_slopes[,3:4]<-d_SMR[,c(5,9)]

      gannotations <- data.frame(
        xpos = c(background_slopes$DateTime_start[1]),
        ypos =  c(max(background_slopes$back_m)+0.5),
        annotateText = c("Used sections of predicted background: \n per L vol"),
        hjustvar = c(0) ,
        vjustvar = c(1.0))

      background_slope_plot<-ggplot(data=back_all, aes(x=DateTime_start, y=m, colour=Ch, group=Ch))+
        geom_point()+
        geom_point(data= background_slopes, mapping = aes(x = DateTime_start, back_m), color = "grey30", size=1, pch=21)+
        geom_line(data= background_slopes, mapping = aes(x = DateTime_start, back_m), color = "grey30", size=1)+
        geom_text(data = gannotations, aes(x=xpos,y=ypos,hjust=hjustvar,
                                           vjust=vjustvar,label=annotateText), color="black")+
        # geom_point(aes(x=DateTime_start, back_m), colour="black", fill="grey", alpha=0.5, pch=19, size=1)+
        facet_grid(Ch~.)+
        # geom_smooth(method="lm", se=FALSE)+
        theme_bw()+
        ylab("Background respiration - mgO2/min/L")+
        theme(axis.text.x = element_text(angle = 45))+
        facet_grid(Ch~.)


		  png(plotname.backgr, width=4, height=8, res=200, units="in")
		    print(background_slope_plot)
		  dev.off()

		  d_SMR$background_slope<- background_slopes$back_m

		}

		if(background_linear_gr==TRUE & match_background_Ch==TRUE){

		  message("SMR/RMR corrected for background: using Ch specific background based on estimated linear regression of bacterial growth")
		  background_slopes<-data.frame(DateTime_start = d_SMR$DateTime_start)

       for (i in 1:nrow(d_SMR)){

		    if(substr(as.character(d_SMR$Ch[i]), start=3, stop=3) == "1"){
		      background_slopes$back_m[i]<-predict(back_regression1, data.frame(DateTime_start = d_SMR$DateTime_start[i]))
          back_m1<-predict(back_regression1, data.frame(DateTime_start = d_SMR$DateTime_start[i]))

  		    d_SMR$mo2[i]<-(( d_SMR$m[i]*(d_SMR$resp.V[i]-d_SMR$bw[i])) - (back_m1 * d_SMR$resp.V[i])) /(d_SMR$bw[i])#^scaling_exponent_smr) # units mgO2 kg-1 min-1 -CORRECTED for back resp
		    # message("Channel: 1,  bacterial respiration slope: ", back_m1)
		    }

		    if(substr(as.character(d_SMR$Ch[i]), start=3, stop=3) == "2"){
		      background_slopes$back_m[i]<-predict(back_regression2, data.frame(DateTime_start = d_SMR$DateTime_start[i]))
		      back_m2<-predict(back_regression2, data.frame(DateTime_start = d_SMR$DateTime_start[i]))
		      d_SMR$mo2[i]<-(( d_SMR$m[i]*(d_SMR$resp.V[i]-d_SMR$bw[i])) - (back_m2 * d_SMR$resp.V[i])) /(d_SMR$bw[i])#^scaling_exponent_smr) # units mgO2 kg-1 min-1 -CORRECTED for back resp
		      # message("Channel: 2,  bacterial respiration slope: ", back_m2)
		    }

		    if(substr(as.character(d_SMR$Ch[i]), start=3, stop=3) == "3"){
		      background_slopes$back_m[i]<-predict(back_regression3, data.frame(DateTime_start = d_SMR$DateTime_start[i]))
		      back_m3<-predict(back_regression3, data.frame(DateTime_start = d_SMR$DateTime_start[i]))
		      d_SMR$mo2[i]<-(( d_SMR$m[i]*(d_SMR$resp.V[i]-d_SMR$bw[i])) - (back_m3 * d_SMR$resp.V[i])) /(d_SMR$bw[i]^scaling_exponent_smr) # units mgO2 kg-1 min-1 -CORRECTED for back resp
		      # message("Channel: 3,  bacterial respiration slope: ", back_m3)
		    }

		    if(substr(as.character(d_SMR$Ch[i]), start=3, stop=3) == "4"){
          background_slopes$back_m[i]<-predict(back_regression4, data.frame(DateTime_start = d_SMR$DateTime_start[i]))
		      back_m4<-predict(back_regression4, data.frame(DateTime_start = d_SMR$DateTime_start[i]))
		      d_SMR$mo2[i]<-(( d_SMR$m[i]*(d_SMR$resp.V[i]-d_SMR$bw[i])) - (back_m4 * d_SMR$resp.V[i])) /(d_SMR$bw[i])#^scaling_exponent_smr) # units mgO2 kg-1 min-1 -CORRECTED for back resp
		      # message("Channel: 4,  bacterial respiration slope: ", back_m4)
		    }
		  }# end of the loop

		  background_slopes[,3:4]<-d_SMR[,c(5,9)]

		  gannotations <- data.frame(
		    xpos = c(background_slopes$DateTime_start[1]),
		    ypos =  c(Inf),
		    annotateText = c("Predicted slopes for back: \n per L vol"),
		    hjustvar = c(0) ,
		    vjustvar = c(1.0))

  		  background_slope_plot<-ggplot(data=back_all, aes(x=DateTime_start, y=m, colour=Ch, group=Ch))+
  		    geom_point()+
  		    geom_point(data= background_slopes, mapping = aes(x = DateTime_start, back_m), color = "grey30", size=1, pch=21)+
  		    geom_line(data= background_slopes, mapping = aes(x = DateTime_start, back_m), color = "grey30", size=1)+
  		    geom_text(data = gannotations, aes(x=xpos,y=ypos,hjust=hjustvar,
  		                                       vjust=vjustvar,label=annotateText), color="black")+
  		    # geom_point(aes(x=DateTime_start, back_m), colour="black", fill="grey", alpha=0.5, pch=19, size=1)+
  		    facet_grid(Ch~.)+
  		    # geom_smooth(method="lm", se=FALSE)+
  		    theme_bw()+
  		    ylab("Regression slope value- mgO2/min/L")+
  		    theme(axis.text.x = element_text(angle = 45))+
  		    facet_grid(Ch~.)

		  png(plotname.backgr, width=4, height=8, res=200, units="in")
		    print(background_slope_plot)
		  dev.off()

		  d_SMR$background_slope<- background_slopes$back_m
		}

		# 1.1 if background files (either prior or post, or both) are provided and its one overall mean value (back_m)
		if ((( !is.null(background_post) | !is.null(background_prior)) & match_background_Ch==FALSE) & is.null(background_slope) & background_linear_gr==FALSE){
		  message("SMR/RMR corrected for background: used a mean (prior and/or post) background measurements | mean bacterial respiration slope: ", back_m, " ~" , round((back_m*100)/(mean(d_SMR[d_SMR$m <= quantile(d_SMR$m, 0.5, na.rm=TRUE), "m"], na.rm=TRUE)), 2), " %")

		  for (i in 1:nrow(d_SMR)){
		    d_SMR$mo2[i]<-(( d_SMR$m[i]*(d_SMR$resp.V[i]-d_SMR$bw[i])) - (back_m * d_SMR$resp.V[i])) /(d_SMR$bw[i])#^scaling_exponent_smr) # units mgO2 kg-1 min-1 -CORRECTED for back resp
		  }

		  d_SMR$background_slope <- back_m
		}

		## if background slope and volume are specifically provided, then use those! this alos overrides the background prior and post argument.
		# all channels with the same slope
		if (!is.null(background_slope)){
		  message("SMR/RMR corrected for background: used a common manually provided background slope for all channels | bacterial respiration slope: ", background_slope, " ~" , round((background_slope*100)/(mean(d_SMR[d_SMR$m <= quantile(d_SMR$m, 0.5, na.rm=TRUE), "m"], na.rm=TRUE)), 2), " %")
		  for (i in 1:nrow(d_SMR)){
		    d_SMR$mo2[i]<-(( d_SMR$m[i]*(d_SMR$resp.V[i]-d_SMR$bw[i])) - (background_slope * background.V)) /(d_SMR$bw[i])#^scaling_exponent_smr) # units mgO2 kg-1 min-1 -CORRECTED for back resp
		  }
		  d_SMR$background_slope<- paste(round(background_slope,2), "_BackVol=", background.V, sep="")
		}

		# 1.2 if background files are provided and its channel specific
		if ((!is.null(background_post) | !is.null(background_prior)) & match_background_Ch==TRUE & background_linear_gr==FALSE){
		  message("SMR/RMR: corrected for background: using Ch specific average background")

       for (i in 1:nrow(d_SMR)){

		    if(substr(d_SMR$Ch[i], start=3, stop=3) == "1"){
  		    d_SMR$mo2[i]<-(( d_SMR$m[i]*(d_SMR$resp.V[i]-d_SMR$bw[i])) - (back_m1 * d_SMR$resp.V[i])) /(d_SMR$bw[i])#^scaling_exponent_smr) # units mgO2 kg-1 min-1 -CORRECTED for back resp
  		    # message("Channel: 1,  bacterial respiration slope: ", back_m1)
  		    d_SMR$background_slope[i]<-back_m1
		    }

		    if(substr(d_SMR$Ch[i], start=3, stop=3) == "2"){
		      d_SMR$mo2[i]<-(( d_SMR$m[i]*(d_SMR$resp.V[i]-d_SMR$bw[i])) - (back_m2 * d_SMR$resp.V[i])) /(d_SMR$bw[i])#^scaling_exponent_smr)# units mgO2 kg-1 min-1 -CORRECTED for back resp
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
		# **********************************************

  	if(any(d_SMR$mo2 < 0)){
      message("!! ATTENTION: negative metabolic rate estimate values")
  	}

  	# size correct each value for SMR calculations, lead with this into EPOC calculations
		d_SMR$mo2_individual_kg<-d_SMR$mo2 * d_SMR$bw
		d_SMR$mo2_common_mass_kg<-exp( (scaling_exponent_smr*(log(d_SMR$common_mass)-log(d_SMR$bw))) + log(d_SMR$mo2_individual_kg))

    # ScalingSlope = (y2 - y1) / (x2 - x1)
    # y2 = log (mo2) <<< current mass dependent value
    # y1 = unknown <<< solve for; this will be log() values! use (exp()) to backtransform
    # x2 = log (mass kg) <<< current body mass in kg
    # x1 = log (common_mass) <<< common body mass for the correction

    # solve for
    # y1 = y2 - ( ScalingSlope * (x2 - x1) )
		# or solve for y2:
		# y2 = ScalingSlope * (x2 - x1) + y1
	  d_SMR$mo2_common_mass_kg <- exp( log(d_SMR$mo2_individual_kg) - (scaling_exponent_smr * (log(d_SMR$bw) - log(d_SMR$common_mass)) ) )
		# b <- exp( (scaling_exponent_smr * ( log(d_SMR$common_mass) - log(d_SMR$bw) ) ) + log(d_SMR$mo2_individual_kg))

	  # rename the mo2 to indicate that this is per kg fish
	  names(d_SMR)[names(d_SMR) == 'mo2'] <- 'mo2_1kg'

    d_SMR$mo2<-d_SMR[,mo2_val_for_calc]

    if(mo2_val_for_calc == "mo2_common_mass_kg"){
       mo2_lab<-bquote(MO[2]~(mgO[2]~min^-1~ .(common_mass) ~ kg^ ~ - .(scaling_exponent_smr)))
       # message(paste("Metabolic performance units: MR/min/ ", common_mass, "kg", sep=""))
      }else{
        if(mo2_val_for_calc == "mo2_individual_kg"){
            mo2_lab<-bquote(MO[2]~(mgO[2]~min^-1~per~individual))
            # message(paste("Metabolic performance units: MR/min/animal"))

          }else{
          	mo2_lab<-bquote(MO[2]~(mgO[2]~min^-1~ kg^-1))
          	# message(paste("Metabolic performance units: MR/min/kg"))
        }
    }

		# find SMR estimates using all methods provided by Chabot
    # --- #
		# 1. frequency plot
		d_SMR$DateTime_start<-as.character(d_SMR$DateTime_start)

		# p_freq<-ggplot(d_SMR, aes(x=mo2))+
		# geom_histogram(bins = 60 , color = "black", fill = "gray") +
		# facet_grid(Ch~.)+
		# theme_classic()+
		# ggtitle(bquote(Datafile: ~ .(data.SMR)))+
		# xlab(mo2_lab)

		# the lowest 10 values after removal of 5 lowest
		a0<-d_SMR[,2:ncol(d_SMR)] %>%
			group_by(Ch)%>%
			top_n(-15, mo2)%>%
			arrange(Ch,mo2)

		# the 5 lowest/ excluded
		a00<-a0 %>%
			group_by(Ch)%>%
			top_n(-5, mo2)%>%
			arrange(Ch,mo2)

		# the 10 lowest/ excluding the 5 lowest in the df
		a<-a0 %>%
			group_by(Ch)%>%
			top_n(10, mo2)%>%
			arrange(Ch,mo2)

		min10_MO2<-as.data.frame(a)

		# detach(package:plyr)
		min10_mean<-min10_MO2%>%
		  group_by(Ch)%>%
			summarize(mean_temp=mean(t_mean), sd_temp=sd(t_mean), mean_mo2=mean(mo2), sd_mo2=sd(mo2),
			cv_mo2 = sd(mo2)/(sqrt(10)), n=length(mo2))

		min10_plot<-ggplot(data=d_SMR, aes(x=min_start, y=mo2))+
			geom_point(size=1)+
			geom_point(data=a00, aes(x=min_start, y=mo2), color="red", pch=19, size=3)+
			geom_point(data=min10_MO2, aes(x=min_start, y=mo2), colour="green4",size=3, alpha=0.7)+
			geom_line(size=0.5, alpha=0.)+
		  theme_light()+
		  ggtitle("Lowest 10 values, excluding the 5 lowest values")+
			theme(legend.position="top")+
			facet_grid(Ch~.)+
		  ylab(mo2_lab)+
		  xlab("Time in trial (min)")


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

		# UPDATE
    colnames(a2)<-c("Ch","quantiles", "mo2_perc")
    a2$mo2_perc<-as.numeric(as.character( a2$mo2_perc))

		quantile_smr<-spread(a2, key=quantiles, value=mo2_perc) # fix
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

		# min15_percPlot<-ggplot(data=d_SMR, aes(x=min_start, y=mo2))+
		# 	geom_point(size=1)+
		# 	geom_point(data=d_SMR[which(d_SMR$Ch=="Ch1" & d_SMR$mo2<=a2$mo2_perc[row_ch1_15perc]),], aes(x=min_start, y=mo2), colour="darkturquoise",size=3)+
		# 	geom_point(data=d_SMR[which(d_SMR$Ch=="Ch2" & d_SMR$mo2<=a2$mo2_perc[row_ch2_15perc]),], aes(x=min_start, y=mo2), colour="darkturquoise",size=3)+
		# 	geom_point(data=d_SMR[which(d_SMR$Ch=="Ch3" & d_SMR$mo2<=a2$mo2_perc[row_ch3_15perc]),], aes(x=min_start, y=mo2), colour="darkturquoise",size=3)+
		# 	geom_point(data=d_SMR[which(d_SMR$Ch=="Ch4" & d_SMR$mo2<=a2$mo2_perc[row_ch4_15perc]),], aes(x=min_start, y=mo2), colour="darkturquoise",size=3)+
		# 	# geom_line(size=0.5, alpha=0.5)+
		# 	theme_classic()+
		# 	ggtitle("15th percentile")+
		# 	theme(legend.position="top")+
		#   ylab(mo2_lab)+
		#   xlab("Time in trial (min)")+
		# 	facet_grid(Ch~.)
		#
		# min20_percPlot<-ggplot(data=d_SMR, aes(x=min_start, y=mo2))+
		# 	geom_point(size=1)+
		# 	geom_point(data=d_SMR[which(d_SMR$Ch=="Ch1" & d_SMR$mo2<=a2$mo2_perc[row_ch1_20perc]),], aes(x=min_start, y=mo2), colour="deepskyblue2",size=3)+
		# 	geom_point(data=d_SMR[which(d_SMR$Ch=="Ch2" & d_SMR$mo2<=a2$mo2_perc[row_ch2_20perc]),], aes(x=min_start, y=mo2), colour="deepskyblue2",size=3)+
		# 	geom_point(data=d_SMR[which(d_SMR$Ch=="Ch3" & d_SMR$mo2<=a2$mo2_perc[row_ch3_20perc]),], aes(x=min_start, y=mo2), colour="deepskyblue2",size=3)+
		# 	geom_point(data=d_SMR[which(d_SMR$Ch=="Ch4" & d_SMR$mo2<=a2$mo2_perc[row_ch4_20perc]),], aes(x=min_start, y=mo2), colour="deepskyblue2",size=3)+
		# 	# geom_line(size=0.5, alpha=0.5)+
		# 	theme_classic()+
		# 	ggtitle("20th percentile")+
		# 	theme(legend.position="top")+
		#   ylab(mo2_lab)+
		#   xlab("Time in trial (min)")+
		# 	facet_grid(Ch~.)
		#
#     if(plot_smr_quantile==10){
# 		  percentile_plot<-min10_percPlot
# 		  }else{
# 		    if(plot_smr_quantile==15){
# 		      percentile_plot<-min15_percPlot
# 		    }else{
# 		      percentile_plot<-min20_percPlot}
# 		}

		png(plotname.freq, width = 10, height = 10, units="in", res=200)
			grid.arrange( min10_plot, min_percPlot, ncol=1, nrow=2)
		dev.off()



		# MLND algorithm from suppl JFB issue 88 Chabot, Steffensen, and Farrell 2016
		# SMR and MMR data frames split based on the channel
		if (length(unique(d_SMR$Ch))>1){
			Ch.data.smr<-split(d_SMR, d_SMR$Ch)
			 	if (!is.null(data.MMR)){
		    	Ch.data.mmr<-split(d_MMR, d_MMR$Ch)
		  	}
		}else{
			Ch.data.smr<-1
		  	if (!is.null(data.MMR)){
		    	Ch.data.mmr<-1
		  	}
		}

		# MLND calculation loop
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
  		if (local_path | !dir.exists("MMR_SMR_AS_EPOC")){
        plotname.mlnd<-paste( filename.SMR,"_", Y.Ch, "_MLND_SMR_analysed.png", sep="")
    	}else{
        plotname.mlnd<-paste("./MMR_SMR_AS_EPOC/plots_mlnd_SMR/", filename.SMR,"_", Y.Ch, "_MLND_SMR_analysed.png", sep="")
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

    			message("SMR/RMR with MLND: only one mean normal distribution to describe data ")
    			if (valid[1]){
    					# message("SMR/RMR with MLND: MLND class 1 = valid")
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
    					# message(paste("SMR/RMR with MLND: lowest valid class = ", valid.clas, sep=""))

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
    					# message(paste("MLND lowest valid class = ", valid.clas, sep=""))

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
    						  # message("More than 3 MLND classes, no diagnostic plots")
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
  	      message("SMR/RMR: MLND method (mean MO2 of the lowest normal distribution) to estimate RMR/SMR measurement is not applied.")
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
  		mmr_overall<-max(mo2)
  		AS_smr_mean10minVal_overall <- mmr_overall - min10_mean$mean_mo2[i]
  		AS_SMR_low10quant_overall <- mmr_overall - quantile_smr[i,2]
  		AS_SMR_low15quant_overall <- mmr_overall - quantile_smr[i,3]
  		AS_SMR_low20quant_overall <- mmr_overall - quantile_smr[i,4]
  		AS_smr_mlnd_overall <- mmr_overall - mlnd

    	values_smr<-as.data.frame(t(c(filename.SMR, ID, Y.Ch, BW, t_min, t_max, t_mean, N_mo2,
  		min10_mean$mean_mo2[i], min10_mean$sd_mo2[i], min10_mean$cv_mo2[i], quantile_smr[i,2],  quantile_smr[i,3], quantile_smr[i,4],
  		mlnd, CVmlnd, Nmlnd,
  		NA, mmr_overall,
  		NA, NA, NA, NA, NA,
  		AS_smr_mean10minVal_overall, AS_SMR_low10quant_overall, AS_SMR_low15quant_overall, AS_SMR_low20quant_overall, AS_smr_mlnd_overall, #5
  		NA, scaling_exponent_mmr, scaling_exponent_smr, common_mass)))#13

  		colnames(values_smr)<-c("filename", "ID", "Ch", "BW","t_min","t_max", "t_mean", "N_mo2", #8
    	"smr_mean10minVal","smr_SD10minVal", "smr_CV10minVal", "SMR_low10quant","SMR_low15quant","SMR_low20quant", #6
    	"smr_mlnd", "smr_CVmlnd", "smr_Nmlnd", #3
    	"mmr", "mmr_overall", #2
    	"AS_smr_mean10minVal", "AS_SMR_low10quant", "AS_SMR_low15quant", "AS_SMR_low20quant", "AS_smr_mlnd", #5
    	"AS_smr_mean10minVal_overall", "AS_SMR_low10quant_overall", "AS_SMR_low15quant_overall", "AS_SMR_low20quant_overall", "AS_smr_mlnd_overall", #5
    	"mmr_length_cycle", "scaling_exponent_mmr", "scaling_exponent_smr", "common_mass")#1
  		newdata.smr<-rbind(newdata.smr, values_smr)

		}# end of the loop for which MLND is a part of



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
		    scale_color_brewer(palette = "Set1")+
				ylab(mo2_lab)+
		    ggtitle("SMR/RMR; # in the circle is (N) mo2 values")+
				xlab("SMR calc type")+
				theme(axis.text.y=element_text(size=20, colour= 'black'),
					axis.text.x=element_text(size=15, colour= 'black', angle=90, hjust=1),
					axis.line.y=element_line(colour = 'black',size=0.5),
					axis.line.x=element_line(colour = 'black',size=0.5),
					axis.ticks.y=element_line(size=0.5),
					axis.ticks.x=element_line(size=0),
					legend.title = element_blank())+
				theme(axis.title.y=element_text(size=15),
					axis.title.x=element_text(size=15),
					panel.border = element_rect(size=0.9,linetype = "solid",fill=NA, colour = "black"))

		png(plotname.smr.meth, width=6, height=5, units="in", res=200)
			print(smr_meth_p)
		dev.off()


	}

	# Both MMR and SMR files are available: recovery focused part
	if(!is.null(data.MMR) & !is.null(data.SMR)){

	  message("Recovery: estimating EPOC, hourly EPOC, percentMMR, time to X percent MMR, and more")

	  if (local_path | !dir.exists("MMR_SMR_AS_EPOC")){
		  	EPOCdata_name<-paste( filename.MMR, "_EPOC_DATA.csv", sep='')
		  	# EPOCplot_name<-	paste( filename.MMR, "_", d$Ch[1], "_EPOC_PLOT.png", sep='')
		  	filename.smr<-paste( gsub('.{4}$', '', data.SMR), "SMR_analyzed.csv", sep='')
  	  	filename.mmr<-paste( gsub('.{4}$', '', data.MMR), "MMR_analyzed.csv", sep='')
	    	filename.MR<-paste( gsub('.{4}$', '', data.MMR), "MR_analyzed.csv", sep='')
		}else{
		  	EPOCdata_name<-paste("./MMR_SMR_AS_EPOC/csv_analyzed_EPOC/", filename.MMR, "_EPOC_DATA.csv", sep='')
		  	# EPOCplot_name<-	paste("../plots_ch_EPOC/", filename.MMR, "_", d$Ch[1], "_EPOC_PLOT.png", sep='')
		  	filename.smr<-paste("./MMR_SMR_AS_EPOC/csv_analyzed_SMR/",gsub('.{4}$', '', data.SMR), "SMR_analyzed.csv", sep='')
  	  	filename.mmr<-paste("./MMR_SMR_AS_EPOC/csv_analyzed_MMR/",gsub('.{4}$', '', data.MMR), "MMR_analyzed.csv", sep='')
	    	filename.MR<-paste("./MMR_SMR_AS_EPOC/csv_analyzed_MR/", gsub('.{4}$', '', data.MMR), "MR_analyzed.csv", sep='')
		}


	  if(!(length(unique(d_SMR$Ch)) == length(unique(d_MMR$Ch)))){
	    stop_function <- TRUE
      if(stop_function) stop(paste("The number of channels in MMR and SMR do not match, check the input data and try again! | "))#, "Ch in MMR:", unique(d_MMR$Ch),  "Ch in SMR:", unique(d_SMR$Ch)))
	  }

		mmr_Ch<-list()
		for(i in 1:4){
			if(any(grepl(as.character(i),as.character(d_MMR$Ch)))){

			  name<-paste("mmr_Ch",i,sep="")
				nameCh<-paste("Ch",i,sep="")
				mmr_Ch0<-d_MMR$mo2[which(d_MMR$cycle_type=="MMR" & d_MMR$Ch==nameCh)]
				mmr_length_cycle0<-d_MMR$cycle_mmr[which(d_MMR$cycle_type=="MMR" & d_MMR$Ch==nameCh)]
				mmr_Ch02<-d_MMR$mo2[which(d_MMR$Ch==nameCh)]
				assign(name, mmr_Ch0)

				mmr_Ch[[i]]<-mmr_Ch0

				cols = c(4:33)
        # newdata.smr[,cols] %<>% lapply(function(x) as.numeric(as.character(x)))
        newdata.smr[,cols] <- lapply(newdata.smr[,cols], as.character)
  	  	newdata.smr[,cols] <- lapply(newdata.smr[,cols], as.numeric)

				newdata.smr$mmr[which(newdata.smr$Ch==nameCh)]<-mmr_Ch0
				newdata.smr$mmr_length_cycle[which(newdata.smr$Ch==nameCh)]<-mmr_length_cycle0
				newdata.smr$mmr_overall[which(newdata.smr$Ch==nameCh)]<-max(c(d_MMR$mo2[which(d_MMR$Ch==nameCh)], d_SMR$mo2[which(d_SMR$Ch==nameCh)]))

				if(mmr_Ch0 < as.numeric(newdata.smr$mmr_overall[which(newdata.smr$Ch==nameCh)])){
				  message(paste("MMR: Spontaneous activity? Fist measurement (", round(mmr_Ch0, 2), ") < highest recorded measurment (", round(max(c(d_MMR$mo2[which(d_MMR$Ch==nameCh)], d_SMR$mo2[which(d_SMR$Ch==nameCh)])), 3)), ")", sep = "")
				}

			}
		}

		# these are lists, so fill in the entire row
		newdata.smr$AS_smr_mean10minVal<-as.numeric(as.character(newdata.smr$mmr))-as.numeric(as.character(newdata.smr$smr_mean10minVal))
		newdata.smr$AS_SMR_low10quant<-as.numeric(as.character(newdata.smr$mmr))-as.numeric(as.character(newdata.smr$SMR_low10quant))
		newdata.smr$AS_SMR_low15quant<-as.numeric(as.character(newdata.smr$mmr))-as.numeric(as.character(newdata.smr$SMR_low15quant))
		newdata.smr$AS_SMR_low20quant<-as.numeric(as.character(newdata.smr$mmr))-as.numeric(as.character(newdata.smr$SMR_low20quant))

		# these are lists, so fill in the entire row
		newdata.smr$AS_smr_mean10minVal_overall<-as.numeric(as.character(newdata.smr$mmr_overall)) -as.numeric(as.character(newdata.smr$smr_mean10minVal))
		newdata.smr$AS_SMR_low10quant_overall<-as.numeric(as.character(newdata.smr$mmr_overall)) -as.numeric(as.character(newdata.smr$SMR_low10quant))
		newdata.smr$AS_SMR_low15quant_overall<-as.numeric(as.character(newdata.smr$mmr_overall)) -as.numeric(as.character(newdata.smr$SMR_low15quant))
		newdata.smr$AS_SMR_low20quant_overall<-as.numeric(as.character(newdata.smr$mmr_overall)) -as.numeric(as.character(newdata.smr$SMR_low20quant))


    if (MLND==TRUE){
      newdata.smr$AS_smr_mlnd<-as.numeric(as.character(newdata.smr$mmr))-as.numeric(as.character(newdata.smr$smr_mlnd))
      newdata.smr$AS_smr_mlnd_overall<-as.numeric(as.character(newdata.smr$mmr_overall)) -as.numeric(as.character(newdata.smr$smr_mlnd))
    }else{
      newdata.smr$AS_smr_mlnd<-0
      newdata.smr$AS_smr_mlnd_overall<-0
    }

    cols = c(4:33)
    cols2 = c(1:3)
  	newdata.smr[,cols] <- lapply(newdata.smr[,cols], as.character)
		newdata.smr[,cols] <- lapply(newdata.smr[,cols], as.numeric)

		newdata.smr[,cols2] <- lapply(newdata.smr[,cols2], as.character)


		# cobine the file specific file to the combined file
		# newdata.all<-rbind(newdata.all, newdata.smr)
  	# EPOC part -- ONLY IF both MMR and SMR are available
  	### defining function first:
  	#	environment(plotEPOC.spar)<-environment()
    # when both MMR and SMR are available

		d_SMR$cycle_type<-"SMR"
		d_MMR$time_frame<-NA

		for(i in 1:nrow(d_MMR)){
			d_MMR$time_frame[i] <- paste("min",d_MMR$cycle_start[i],"_", d_MMR$cycle_end[i], sep="")
		}

		names(d_MMR)[names(d_MMR) == 'cycle_start'] <- 'min_start'

		dat_MMR<-d_MMR[,c("ID", "time_frame", "min_start", "r2", "b", "m", "t_min", "t_max", "t_mean", "Ch", "bw", "mo2", "cycle_type", "DateTime_start", "scaling_exponent", "common_mass")]
		dat_SMR<-d_SMR[,c("ID", "time_frame", "min_start", "r2", "b", "m", "t_min", "t_max", "t_mean", "Ch", "bw", "mo2", "cycle_type", "DateTime_start", "scaling_exponent", "common_mass")]

		if (length(unique(dat_SMR$Ch))>1){
			Ch.dat.smr<-split(dat_SMR, dat_SMR$Ch)
			Ch.dat.mmr<-split(dat_MMR, dat_MMR$Ch)
		}else{
			Ch.dat.smr<-1
			Ch.dat.mmr<-1
		}
#~ 		data<-dat_MMR[-c(1:nrow(dat_MMR)), ]

		EPOCdata<-matrix(ncol=26, nrow=0)
		colnames(EPOCdata) <- c("ID", "smr_type", "smr","spar", "EPOC_full", "end_EPOC_min", "SMR_intergral_full", "SMR_threshold", "EPOC_1hr", "MO2_1hr", "EPOC_2hr", "MO2_2hr", "EPOC_3hr", "MO2_3hr", "EPOC_4hr", "MO2_4hr",  "EPOC_5hr", "MO2_5hr", "end_EPOC.mmr", "EPOC_mmr", "MO2_mmr", "MMR", "MMR_percent", "scaling_exponent_mmr",  "scaling_exponent_smr", "common_mass")

    # apply EPOC.spar function on combined mmr smr data
  	# loop through all available chanels
		for(i in 1:length(unique(dat_SMR$Ch))){

  		if (length(unique(dat_SMR$Ch))==1){
  			d.smr<-dat_SMR
  			d.mmr<-dat_MMR
  		}else{
  			d.smr<-as.data.frame(Ch.dat.smr[i])
  			d.mmr<-as.data.frame(Ch.dat.mmr[i])
  		}


		  cols = c(1,2,5,13,14)
		  d.smr[,cols] <- lapply(d.smr[, cols], as.character)
      d.mmr[,cols] <- lapply(d.mmr[, cols], as.character)

  		d<-rbind(d.mmr, d.smr)
  		colnames(d)<-c("ID", "time_frame", "min_start", "r2", "b", "m", "t_min", "t_max", "t_mean", "Ch", "bw", "mo2", "cycle_type", "DateTime_start", "scaling_exponent", "common_mass")

  		# mmr value for MMR 50 recovery calculations
  		mmr.val <- max(d$mo2)
  		# message(paste(d$Ch[1], ": MMR = ", round(mmr.val,2) ,sep=""))

  		# Ch specific name for epoc plot
			if (local_path | !dir.exists("MMR_SMR_AS_EPOC")){
	    	EPOCplot_name<-	paste( filename.MMR, "_", d$Ch[1], "_EPOC_PLOT.png", sep='')
    	}else{
	    	EPOCplot_name<-	paste("./MMR_SMR_AS_EPOC/plots_ch_EPOC/", filename.MMR, "_", d$Ch[1], "_EPOC_PLOT.png", sep='')
	    }

    	 # d$DateTime_start<- strptime(d$DateTime_start, format="%Y-%m-%d %H:%M:%S") # the default strptime, that was already used above
  		 d$time_mo2<-NA
  		 d$time_mo2[1]<-0

			for(j in 2:nrow(d)){
				d$time_mo2[j]<-difftime(d$DateTime_start[j], d$DateTime_start[1], units=c("mins"))
			}

		  end_EPOC<-end_EPOC_Ch[as.numeric(substr(d$Ch[1], start=3, stop=3))]
  		## The EPOC calculation

		  spars <- c(0.1, 0.2, 0.3)

  			if(nrow(d) ==4 || nrow(d)>4){
  				for (n in 1:length(spars)){
  					if (n == 1) {
  					  # length(spars)
  						png(EPOCplot_name, width = 6, height = length(spars)*4, units="in", res=200)
  						par(mfrow=c(length(spars),1), mar=c(4,4,3,1)+0.1)
  					}

  					EPOCdata <- EPOC.spar(spars[n], d, EPOCdata, mmr.val, epoc_threshold, recovMMR_threshold, newdata.smr, MLND, end_EPOC, scaling_exponent_mmr, scaling_exponent_smr, common_mass)

  					if (n == length(spars)){
  						dev.off()
  						write.csv(file=EPOCdata_name, EPOCdata, row.names=FALSE)
  					}
  				}

  			}else{
  				message(paste("Recovery: Too few (n=", nrow(d), ") data points for smoothing and EPOC analysis.", sep = ""))
  			}

	   }


  	if (epoc_threshold == 1){
  		  	message (paste("Recovery: Time at EPOC assigned when oxygen uptake rates recover to RMR/SMR level"))
  		}else{
  		    message (paste("Recovery: Time at EPOC assigned when oxygen uptake rates recover to ", epoc_threshold, " x SMR/RMR ", sep = ""))
  		}

#   	newdata.smr<-as.data.frame(newdata.smr)
# 		newdata.smr<- apply(newdata.smr, 2, as.character)
#
		write.csv(file=filename.smr, d_SMR, row.names=FALSE)
		write.csv(file=filename.mmr, d_MMR, row.names=FALSE)
		write.csv(file=filename.MR, newdata.smr, row.names=FALSE)

	# 	if(test=="preSDA"){
	#     	  filename.MR<-paste("../csv_input_files/", gsub('.{4}$', '', data.MMR), "MR_analyzed.csv", sep='')
	#     	  write.csv(file=filename.MR, newdata.smr, row.names=FALSE)
	#    }

		# message("Save MMR, SMR, and MR files")

  }

	# no MMR file, only SMR file
  if(!is.null(data.SMR) & is.null(data.MMR)){

   message("No RMR/SMR: Aerobic scopes and recovery performances are not estimated.")

	 if (local_path | !dir.exists("MMR_SMR_AS_EPOC")){
	    filename.smr<-paste(gsub('.{4}$', '', data.SMR), "SMR_analyzed.csv", sep='')
		  filename.MR<-paste(gsub('.{4}$', '', data.SMR), "MR_analyzed.csv", sep='')
   }else{
	    filename.smr<-paste("./MMR_SMR_AS_EPOC/csv_analyzed_SMR/",gsub('.{4}$', '', data.SMR), "SMR_analyzed.csv", sep='')
		  filename.MR<-paste("./MMR_SMR_AS_EPOC/csv_analyzed_MR/", gsub('.{4}$', '', data.SMR), "MR_analyzed.csv", sep='')
	 }

    lst <- lapply(newdata.smr, unlist)
    newdata.smr <- (data.frame(lapply(lst, `length<-`, max(lengths(lst)))))

		write.csv(file=filename.smr, d_SMR, row.names=FALSE)
		write.csv(file=filename.MR, newdata.smr, row.names=FALSE)

	   # if(test=="preSDA"){
	   #  	  filename.MR<-paste("../csv_input_files/", gsub('.{4}$', '', data.SMR), "MR_analyzed.csv", sep='')
	   #  	  write.csv(file=filename.MR, newdata.smr, row.names=FALSE)
	   # }
		# message("Save SMR and MR files")

  }


	# ***********************************************
  # no SMR file
  if(is.null(data.SMR) & !is.null(data.MMR)){# if only MMR is analyzed

    message("No MMR: Aerobic scopes and recovery performances are not estimated.")

		# filename.MMR<-paste(gsub('.{3}$', '',data.MMR), "_MMR_", sep="")
	  if (local_path | !dir.exists("MMR_SMR_AS_EPOC")){
	    filename.mmr<-paste( gsub('.{4}$', '', data.MMR), "MMR_analyzed.csv", sep='')
		  filename.MR<-paste( gsub('.{4}$', '', data.MMR), "MR_analyzed.csv", sep='')

	   }else{
	    filename.mmr<-paste("./MMR_SMR_AS_EPOC/csv_analyzed_MMR/",gsub('.{4}$', '', data.MMR), "MMR_analyzed.csv", sep='')
		  filename.MR<-paste("./MMR_SMR_AS_EPOC/csv_analyzed_MR/", gsub('.{4}$', '', data.MMR), "MR_analyzed.csv", sep='')
    }

		# d_MMR$bw<-NA
		# d_MMR$mo2<-NA
		# d_MMR$ID<-NA
		# d_MMR$resp.V<-NA
		#
		# for(i in 1:4){
		# 	if(any(grepl(as.character(i),as.character(d_MMR$Ch)))){
		# 		nameCh<-paste("Ch",i,sep="")
		# 		bw.val<-BW.animal[i]
		# 		ID<-AnimalID[i]
		# 		resp.Vol<-resp.V[i]
		# 		n.row<-which(d_MMR$Ch==nameCh)
		# 		d_MMR$bw[n.row]<-bw.val
		# 		d_MMR$ID[n.row]<-ID
		#     d_MMR$resp.V[n.row]<-resp.Vol
		# 	}
		# }
		#
	# ### CORERCT FOR background!!!
	# 	for (i in 1:nrow(d_MMR)){
	# 		d_MMR$mo2[i]<-d_MMR$m[i]*(d_MMR$resp.V[i]-d_MMR$bw[i])/(d_MMR$bw[i])#^scaling_exponent_mmr) # units mgO2 kg-1 min-1
	# 	}

		d_MMR<-d_MMR[as.character(d_MMR$ID)!="blank",]
	  d_MMR<-d_MMR[d_MMR$ID!="NA",]


		for (i in 1:nrow(d_MMR[d_MMR$cycle_type=="MMR",])){
		  d_MMR2<- d_MMR[d_MMR$cycle_type=="MMR",]# data frame without the cycles (e.g., the cycles only

		  # temp min and max for MMR only
			values<-as.data.frame(t(c(gsub('.{4}$', '', data.MMR),  d_MMR2$ID[i], d_MMR2$Ch[i],
			                          d_MMR2$bw[i], d_MMR2$t_min[i], d_MMR2$t_max[i], d_MMR2$t_mean[i], NA,
			NA,NA,NA,
			NA,NA,NA,
			NA,NA,NA,
			d_MMR2$mo2[i],NA,NA,NA,
			NA,NA,NA,
			NA,NA,NA,
			NA,NA, d_MMR2$cycle_mmr[i],
			scaling_exponent_mmr, NA, d_MMR2$common_mass[i])))

			colnames(values)<-c("filename", "ID", "Ch",
			                    "BW","t_min","t_max", "t_mean", "N_mo2",
			"smr_mean10minVal","smr_SD10minVal", "smr_CV10minVal",
			"SMR_low10quant","SMR_low15quant",  "SMR_low20quant",
			"smr_mlnd", "smr_CVmlnd", "smr_Nmlnd",
			"mmr" ,"mmr_overall" , "AS_smr_mean10minVal", "AS_SMR_low10quant",
			"AS_SMR_low15quant" , "AS_SMR_low20quant", "AS_smr_mlnd",
			"AS_smr_mean10minVal_overall","AS_SMR_low10quant_overall", "AS_SMR_low15quant_overall",
			"AS_SMR_low20quant_overall","AS_smr_mlnd_overall", "mmr_length_cycle",
			"scaling_exponent_mmr", "scaling_exponent_smr", "common_mass")

			newdata.smr<-rbind(newdata.smr, values)
		}




		write.csv(file=filename.mmr, d_MMR, row.names=FALSE)

		# newdata.smr<-as.data.frame(newdata.smr)
		# newdata.smr<- apply(newdata.smr, 2, as.character)
		lst <- lapply(newdata.smr, unlist)
    newdata.smr <- (data.frame(lapply(lst, `length<-`, max(lengths(lst)))))

		write.csv(file=filename.MR, newdata.smr, row.names=FALSE)
	# 	if(test=="preSDA"){
	#     	  filename.MR<-paste("../csv_input_files/", gsub('.{4}$', '', data.MMR), "MR_analyzed.csv", sep='')
	#     	  write.csv(file=filename.MR, newdata.smr, row.names=FALSE)
	#   }
		# message("Saving MMR, and MR files")
  }

}









