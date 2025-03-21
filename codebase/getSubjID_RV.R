#-------------------------------------------------------
# Author: Brady M. Chisholm
# University of Minnesota Twin Cities, Dpt. of Psychology
# 3.17.2025
#-------------------------------------------------------

# NOTE: subjID is not subject order, nor is it chronological. Utilize: 
# ln 10 = subj 1-10; ln 11 = 11-20; ln 12 = subj 16-30; ln 13 = subj 31-46; ln 14 = subj 47-

getSubjList <- function(experiment) {
if (experiment == 1) {
  subjList <- c('P01', 'P02', 'P03', 'P04', 'P05', 'P06', 'P07', 'P08', 'P09', 'P10',
                'P11', 'P12', 'P13', 'HI14', 'HI15', 'P17', 'HI18', 'P19', 'HI20', 'HI21',
                'P22', 'HI23', 'P24', 'P25', 'P26', 'HI27', 'HI28', 'P30', 'P31', 'P32',
                'HI33', 'HI34', 'P35', 'HI36', 'P37', 'P38', 'HI39', 'P40', 'HI41', 'HI42',
                'HI43', 'HI44', 'HI45', 'HI46', 'HI47', 'HI48', 'HI51', 'P41', 'HI50', 'HI52')
  
excludeTheseSubjs <- c('P01','P02', 'P35', 'HI49')   # some bug with P02 subjID data, exclude for now 

# P01 incomplete; P35 no main session; HI49 has detached retina.
# Subj XX excluded in original full data b/c of EEG, but o.k. pupil data.
subjList <- subjList[!(subjList %in% excludeTheseSubjs)]
} else {
  # other exp subj data `placeholder`
  print("ERROR: Input given for which data does not exist.")
  subjList <- NULL
}
return(subjList)
}


