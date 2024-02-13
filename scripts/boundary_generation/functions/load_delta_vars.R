library(readxl)


# load_delta_vars: Function -----------------------------------------------
# reads in a variable from the DSM2 dataset


load_delta_vars <- function(dsm2_filename, var_out) {
  sheetlist <- excel_sheets(path=dsm2_filename)
  
  for (sheet in sheetlist) {
    dsm2_inouts <- as.data.frame(read_excel(dsm2_filename, sheet=sheet))
    names(dsm2_inouts)[1] <- 'Time'
    dsm2_inouts$Time <- lubridate::ymd(dsm2_inouts$Time)
    
    if (sheet==var_out) {
      return(dsm2_inouts)
    }
  }
}