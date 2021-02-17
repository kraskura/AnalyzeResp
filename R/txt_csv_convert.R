
#' Title
#'
#' @param txt_file the raw text file from Firesting
#' @param N_Ch Number of channels
#' @param path the paths where to save the csv file
#'
#' @return generated csv file, automatically saved

# txt_csv_convert<-function(txt_file, N_Ch = 4, path = ".", save = FALSE){
#
#   if(N_Ch == 4 | N_Ch==2){
#     new_csv<-as.data.frame(matrix(nrow=0, ncol=8))
#     colnames(new_csv)<-c("date", "time", "time_sec", "Ch1_O2", "Ch1_temp", "Ch2_O2", "Ch3_O2", "Ch4_O2")
#
#     d<-read.delim(txt_file, skip=19)
#
#
#     # don't think need this anymore
#     if (!colnames(d[1])=="Date"){
#       d<-read.delim(txt_file, skip=26)
#     }
#
#     nr<-nrow(d)
#     nc<-ncol(d)
#
#     new_csv[nr,]<-NA
#     new_csv$date<-d[,1]
#     new_csv$time<-d[,2]
#     new_csv$time_sec<-d[,3]
#     # oxyegen below
#     new_csv$Ch1_O2<-d[,5]
#
#     new_csv$Ch1_temp<-d[,15]# temp Ch1 - but same for all
#
#     new_csv$Ch2_O2<-d[,6]
#     new_csv$Ch3_O2<-d[,7]
#     new_csv$Ch4_O2<-d[,8]
#   } # end od N_Ch == 2
#
#
#   if (N_Ch==8){
#     new_csv<-as.data.frame(matrix(nrow=0, ncol=11))
#     colnames(new_csv)<-c("date", "time", "time_sec", "Ch1_O2", "Ch1_temp", "Ch2_O2", "Ch3_O2", "Ch4_O2", "Ch2_temp", "Ch3_temp", "Ch4_temp")
#
#     d<-read.delim(txt_file, skip=26)
#
#     nr<-nrow(d)
#     nc<-ncol(d)
#
#     new_csv[nr,]<-NA
#     new_csv$date<-d[,1]
#     new_csv$time<-d[,2]
#     new_csv$time_sec<-d[,3]
#     # oxyegen below
#     new_csv$Ch1_O2<-d[,5]
#
#     new_csv$Ch1_temp<-d[,9]# unique ch temps
#     new_csv$Ch2_temp<-d[,10]# unique ch temps
#     new_csv$Ch3_temp<-d[,11]# unique ch temps
#     new_csv$Ch4_temp<-d[,12]# unique ch temps
#
#     new_csv$Ch2_O2<-d[,6]
#     new_csv$Ch3_O2<-d[,7]
#     new_csv$Ch4_O2<-d[,8]
#   }
#
#
#   if(save==TRUE){
#     # save in the current directory (default)
#     if(path == "."){
#       write.csv(file=paste(gsub('.{4}$', '', txt_file), ".csv", sep=''), new_csv, row.names=FALSE)
#     }
#
#     #save in newly created folder - AUTO; use function: organize_MR_analysis()
#     if(path == "AUTO"){
#       write.csv(file=paste("AUTO/csv_files/", gsub('.{4}$', '', txt_file), ".csv", sep=''), new_csv, row.names=FALSE)
#     }
#
#     #save in newly created folder - MANUAL; organize_MR_analysis()
#     if(path == "MANUAL"){
#       write.csv(file=paste("MANUAL/csv_files/", gsub('.{4}$', '', txt_file), ".csv", sep=''), new_csv, row.names=FALSE)
#     }
#
#     # save in newly created folder - BACTERIAL_RESP; organize_MR_analysis()
#     if(path == "BACTERIAL_RESP"){
#       write.csv(file=paste("BACTERIAL_RESP/csv_files/", gsub('.{4}$', '', txt_file), ".csv", sep=''), new_csv, row.names=FALSE)
#     }
#
#     if(path == "SDA"){
#       write.csv(file=paste("SDA/csv_files/", gsub('.{4}$', '', txt_file), ".csv", sep=''), new_csv, row.names=FALSE)
#     }
#   }
#
# }
