#' Title
#'
#' @return The output from \code{\link{print}}
#' @export
organizeAnalysisLocally <- function (){

  # userConfirm<-readline(prompt="This function will create a set of data folder at a local directory: (y/n) ")
  # if(userConfirm == "Y" | userConfirm == "y" | userConfirm == "n" | userConfirm == "N"){

    # if(userConfirm == "n" | userConfirm == "N"){
      # message("Folders not created, all files need to be handled within local directory (default for all functions)")
      # }else{
      userInputFolders<-readline(prompt="Indicate all folders to create e.g. 1,2,3: \n (1) SMR \n (2) MMR \n (3) SDA \n (4) MMR_SMR_AS_EPOC \n (5) BACK_RESP")
      userInputFolders<-(c(as.numeric(unlist(strsplit(userInputFolders,split=',')))))

      if (!dir.exists("csv_files")){
        dir.create(file.path("./csv_files"), recursive = TRUE)
      }

      if(any(userInputFolders == 1)){
        if (!dir.exists("./SMR")){
          dir.create(file.path("SMR"), recursive = TRUE)
          dir.create(file.path("SMR","csv_analyzed"), recursive = TRUE)
          dir.create(file.path("SMR","plots_all_together"), recursive = TRUE)
          dir.create(file.path("SMR","plots_channel"), recursive = TRUE)
          dir.create(file.path("SMR","plots_summary_respo"), recursive = TRUE)
        }
      }

      if(any(userInputFolders == 2)){
        if (!dir.exists("./MMR")){
          dir.create(file.path("MMR"), recursive = TRUE)
          dir.create(file.path("MMR","csv_analyzed"), recursive = TRUE)
          dir.create(file.path("MMR","channel_plots"), recursive = TRUE)
          dir.create(file.path("MMR","channel_plots_MMRanalysis"), recursive = TRUE)
          dir.create(file.path("MMR","channel_sliding_sets"), recursive = TRUE)
        }
      }

      if(any(userInputFolders == 3)){
         if (!dir.exists("./SDA")){
          dir.create(file.path("SDA"), recursive = TRUE)
          dir.create(file.path("SDA", "csv_analyzed_SDA"), recursive = TRUE) # equivalent to EPOC for the MMR_SMR_AS_EPOC function
          dir.create(file.path("SDA", "csv_analyzed_SMR"), recursive = TRUE) # equivalent to EPOC for the MMR_SMR_AS_EPOC function
          dir.create(file.path("SDA", "csv_analyzed_SDA_hrly"), recursive = TRUE)
          dir.create(file.path("SDA", "csv_analyzed_MR"), recursive = TRUE)
          dir.create(file.path("SDA", "csv_input_files"), recursive = TRUE)


          dir.create(file.path("SDA", "csv_analyzed"), recursive = TRUE)
          dir.create(file.path("SDA", "plots_ch_SDA"), recursive = TRUE) # equivalent to EPOC for the MMR_SMR_AS_EPOC function
          dir.create(file.path("SDA", "plots_methods_sum_SMR"), recursive = TRUE)
          dir.create(file.path("SDA", "plots_min_values_SMR"), recursive = TRUE)
          dir.create(file.path("SDA", "plots_mlnd_SMR"), recursive = TRUE)
          dir.create(file.path("SDA", "plots_SDA_hourly"), recursive = TRUE)
          dir.create(file.path("SDA","plots_summary_respo"), recursive = TRUE)
          dir.create(file.path("SDA","plots_channel"), recursive = TRUE)
          dir.create(file.path("SDA","plots_all_together"), recursive = TRUE)
         }
      }

      if(any(userInputFolders == 4)){

       if (!dir.exists("./MMR_SMR_AS_EPOC")){
          dir.create(file.path("MMR_SMR_AS_EPOC"), recursive = TRUE)
          dir.create(file.path("MMR_SMR_AS_EPOC", "csv_analyzed_EPOC"), recursive = TRUE)
          dir.create(file.path("MMR_SMR_AS_EPOC", "csv_analyzed_MMR"), recursive = TRUE)
          dir.create(file.path("MMR_SMR_AS_EPOC", "csv_analyzed_SMR"), recursive = TRUE)
          dir.create(file.path("MMR_SMR_AS_EPOC", "csv_analyzed_MR"), recursive = TRUE)
          dir.create(file.path("MMR_SMR_AS_EPOC", "csv_input_files"), recursive = TRUE)
          dir.create(file.path("MMR_SMR_AS_EPOC", "plots_ch_EPOC"), recursive = TRUE)
          dir.create(file.path("MMR_SMR_AS_EPOC", "plots_methods_sum_SMR"), recursive = TRUE)
          dir.create(file.path("MMR_SMR_AS_EPOC", "plots_min_values_SMR"), recursive = TRUE)
          dir.create(file.path("MMR_SMR_AS_EPOC", "plots_mlnd_SMR"), recursive = TRUE)
       }
      }

      if(any(userInputFolders == 5)){
        if (!dir.exists("./BACK_RESP")){
          dir.create(file.path("BACK_RESP"), recursive = TRUE)
          dir.create(file.path("BACK_RESP","csv_analyzed"), recursive = TRUE)
          dir.create(file.path("BACK_RESP","plots_all_together"), recursive = TRUE)
          dir.create(file.path("BACK_RESP","plots_channel"), recursive = TRUE)
          dir.create(file.path("BACK_RESP","plots_summary_respo"), recursive = TRUE)
        }
      }
    # }
  # }else{
  #     message("Unexpected input: Folders not created, all files need to be handled within local directory (default for all functions)")
  # }
}
