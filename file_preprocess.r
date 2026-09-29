# function for data rearrangement
# aggregations and logits
#
# library(this.path)
# setwd(this.dir())


# read file and duplicate and sort columns for easy rearrangement
file_preprocess <- function(datafile = "../Data/NOVICE_rawDataAnonymized.xlsx"){
  library(readxl)
  daten <- read_excel(datafile)
  k <- ncol(daten)
  daten$Autonomie2 <- daten$Autonomie
  daten$Autonomie3 <- daten$Autonomie
  daten$AutonomieE2 <- daten$AutonomieE
  daten$AutonomieE3 <- daten$AutonomieE
  daten$AutonomieD2 <- daten$AutonomieD
  daten$AutonomieD3 <- daten$AutonomieD
  daten <- daten[, c(1, 4:9,                 # case_vars   1:7
                     56:68, k+5, k+6, 69:71, # case_sums   8:25
                     2, 10:16,               # rep_vars1  26:33
                     30:32,                  # rep_sums1  34:36
                     17:29, k+1, k+2,        # item_vars1 37:51
                     3, 33:39,               # rep_vars2  52:59
                     53:55,                  # rep_sums2  60:62
                     40:52, k+3, k+4)]       # item_vars2 63:77
  gesamt <- substr(names(daten), 1, 6) == "Gesamt"
  names(daten)[gesamt] <- paste0("GOALS", substr(names(daten)[gesamt], 7, 8))
  return(as.data.frame(daten))
} # file_preprocess
daten <- file_preprocess(datafile = "../Data/NOVICE_rawDataAnonymized.xlsx")
# names(daten)

