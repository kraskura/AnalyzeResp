#' @title convenient function to create a set of local working directories to streamline workflow
#'
#' @description
#' Creates a collection of local directories that help organizing the analysis by exporting output files to specific folders associated with each metabolic performance function (not required for smooth data analysis).
#'
#' @param SMR.folder create a local SMR folder, default = TRUE
#' @param MMR.folder create a local MMR folder, default = TRUE
#' @param SDA.folder create a local SDA folder, default = TRUE
#' @param MMR_SMR_AS_EPOC.folder create a local MMR_SMR_AS_EPOC.folder, default = TRUE
#' @param BACK_RESP.folder create a local BACK_RESP.folder, default = TRUE
#'
#' @return The output from \code{\link{print}}
#' @export
organizeAnalysisLocally <- function (SMR.folder = TRUE,
                                     MMR.folder = TRUE,
                                     SDA.folder = TRUE,
                                     MMR_SMR_AS_EPOC.folder = TRUE,
                                     BACK_RESP.folder = TRUE){

     message("Package was last updated: ",  Sys.time())
      # userInputFolders<-readline(prompt="Indicate all folders to create e.g. type: 1, 2, 3: (1) SMR (2) MMR (3) SDA (4) MMR_SMR_AS_EPOC (5) BACK_RESP")
      # userInputFolders<-(c(as.numeric(unlist(strsplit(userInputFolders,split=',')))))

      if (!dir.exists("csv_files")){
        dir.create(file.path("./csv_files"), recursive = TRUE)
      }

      if(SMR.folder){
        if (!dir.exists("./SMR")){
          dir.create(file.path("SMR"), recursive = TRUE)
          dir.create(file.path("SMR","csv_analyzed"), recursive = TRUE)
          dir.create(file.path("SMR","plots_all_together"), recursive = TRUE)
          dir.create(file.path("SMR","plots_channel"), recursive = TRUE)
          dir.create(file.path("SMR","plots_summary_respo"), recursive = TRUE)
        }
      }

      if(MMR.folder){
        if (!dir.exists("./MMR")){
          dir.create(file.path("MMR"), recursive = TRUE)
          dir.create(file.path("MMR","csv_analyzed"), recursive = TRUE)
          dir.create(file.path("MMR","channel_plots"), recursive = TRUE)
          dir.create(file.path("MMR","channel_plots_MMRanalysis"), recursive = TRUE)
          dir.create(file.path("MMR","channel_sliding_sets"), recursive = TRUE)
        }
      }

      if(SDA.folder){
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

      if(MMR_SMR_AS_EPOC.folder){

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

      if(BACK_RESP.folder){
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
