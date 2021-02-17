
#' Title
#'
#' @param create folders to create to organise local export datafiles
#'
#' @return None, saved folder
#' @export
#'
#' @examples
#' \dontrun{
#' organize_MR_analysis(create="none")
#' }
organize_MR_analysis <- function (create = c("AUTO", "MANUAL", "Full", "BACTERIAL_RESP", "SDA")) {

  if(create == "AUTO"){
    if (!dir.exists("./AUTO")){
      dir.create(file.path("AUTO"), recursive = TRUE)
      dir.create(file.path("AUTO","csv_analyzed"), recursive = TRUE)
      dir.create(file.path("AUTO","csv_files"), recursive = TRUE)
      dir.create(file.path("AUTO","plots_all_together"), recursive = TRUE)
      # dir.create(file.path("AUTO","plots_channel_temperature"), recursive = TRUE)
      dir.create(file.path("AUTO","plots_channel_cycle"), recursive = TRUE)
      dir.create(file.path("AUTO","plots_channel"), recursive = TRUE)
      dir.create(file.path("AUTO","plots_summary_respo"), recursive = TRUE)
    }
  }

  if(create == "MANUAL"){
    if (!dir.exists("./MANUAL")){
      dir.create(file.path("MANUAL"), recursive = TRUE)
      dir.create(file.path("MANUAL","csv_analyzed"), recursive = TRUE)
      dir.create(file.path("MANUAL","csv_files"), recursive = TRUE)
      dir.create(file.path("MANUAL","channel_plots"), recursive = TRUE)
      dir.create(file.path("MANUAL","channel_plots_MMRanalysis"), recursive = TRUE)
      dir.create(file.path("MANUAL","channel_sliding_sets"), recursive = TRUE)
    }
  }

  if(create == "BACTERIAL_RESP"){
    if (!dir.exists("./BACTERIAL_RESP")){
      dir.create(file.path("BACTERIAL_RESP"), recursive = TRUE)
      dir.create(file.path("BACTERIAL_RESP","csv_analyzed"), recursive = TRUE)
      dir.create(file.path("BACTERIAL_RESP","csv_files"), recursive = TRUE)
      dir.create(file.path("BACTERIAL_RESP","plots_all_together"), recursive = TRUE)
      # dir.create(file.path("BACTERIAL_RESP","plots_channel_temperature"), recursive = TRUE)
      dir.create(file.path("BACTERIAL_RESP","plots_channel_cycle"), recursive = TRUE)
      dir.create(file.path("BACTERIAL_RESP","plots_channel"), recursive = TRUE)
      dir.create(file.path("BACTERIAL_RESP","plots_summary_respo"), recursive = TRUE)
    }
  }

  if(create == "SDA"){
    if (!dir.exists("./SDA")){
      dir.create(file.path("SDA"), recursive = TRUE)
      dir.create(file.path("SDA", "csv_analyzed_SDA"), recursive = TRUE) # equivalent to EPOC for the MMR_SMR_AS_EPOC function
      dir.create(file.path("SDA", "csv_analyzed_SMR"), recursive = TRUE) # equivalent to EPOC for the MMR_SMR_AS_EPOC function
      dir.create(file.path("SDA", "csv_analyzed_SDA_hrly"), recursive = TRUE)
      dir.create(file.path("SDA", "csv_analyzed_MR"), recursive = TRUE)
      dir.create(file.path("SDA", "csv_input_files"), recursive = TRUE)


      dir.create(file.path("SDA", "csv_analyzed"), recursive = TRUE)
      dir.create(file.path("SDA", "csv_files"), recursive = TRUE)
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

  if(create == "Full"){

    if (!dir.exists("./AUTO")){
      dir.create(file.path("AUTO"), recursive = TRUE)
      dir.create(file.path("AUTO","csv_analyzed"), recursive = TRUE)
      dir.create(file.path("AUTO","csv_files"), recursive = TRUE)
      dir.create(file.path("AUTO","plots_all_together"), recursive = TRUE)
      # dir.create(file.path("AUTO","plots_channel_temperature"), recursive = TRUE)
      dir.create(file.path("AUTO","plots_channel_cycle"), recursive = TRUE)
      dir.create(file.path("AUTO","plots_channel"), recursive = TRUE)
      dir.create(file.path("AUTO","plots_summary_respo"), recursive = TRUE)
    }

    if (!dir.exists("./BACTERIAL_RESP")){
      dir.create(file.path("BACTERIAL_RESP"), recursive = TRUE)
      dir.create(file.path("BACTERIAL_RESP","csv_analyzed"), recursive = TRUE)
      dir.create(file.path("BACTERIAL_RESP","csv_files"), recursive = TRUE)
      dir.create(file.path("BACTERIAL_RESP","plots_all_together"), recursive = TRUE)
      # dir.create(file.path("BACTERIAL_RESP","plots_channel_temperature"), recursive = TRUE)
      dir.create(file.path("BACTERIAL_RESP","plots_channel_cycle"), recursive = TRUE)
      dir.create(file.path("BACTERIAL_RESP","plots_channel"), recursive = TRUE)
      dir.create(file.path("BACTERIAL_RESP","plots_summary_respo"), recursive = TRUE)
    }

    if (!dir.exists("./MANUAL")){
      dir.create(file.path("MANUAL"), recursive = TRUE)
      dir.create(file.path("MANUAL","csv_analyzed"), recursive = TRUE)
      dir.create(file.path("MANUAL","csv_files"), recursive = TRUE)
      dir.create(file.path("MANUAL","channel_plots"), recursive = TRUE)
      dir.create(file.path("MANUAL","channel_plots_MMRanalysis"), recursive = TRUE)
      dir.create(file.path("MANUAL","channel_sliding_sets"), recursive = TRUE)
    }

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
}
