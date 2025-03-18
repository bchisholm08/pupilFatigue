#-------------------------------------------------------
# Author: Brady M. Chisholm
# University of Minnesota Twin Cities, Dpt. of Psychology
# Date: 3.17.2025
#
# Description: This script loads raw MATLAB pupil data, checks and (if necessary)
#              upsamples the data to 300Hz (6,000 samples for 20 sec trials),
#              downsampled to 150Hz (3,000 samples), selects the better eye,
#              and saves the processed data.
#-------------------------------------------------------

# --- Integrated Resample Signal Function ---
if (!require(signal)) {
  install.packages("signal")
  library(signal)
}

# Debugged version of resample_signal function
resample_signal <- function(inputMatrix, upSamp, downSamp, filt = 10, beta = 5) { 
  # If inputMatrix is a vector, convert to a one-column matrix for uniform handling.
  original_dim <- dim(inputMatrix)
  if (is.null(original_dim)) {
    inputMatrix <- matrix(inputMatrix, ncol = 1)
  }
  
  # Determine filter length:
  factor <- max(1, downSamp / upSamp)
  L <- as.integer(2 * filt * factor + 1)
  
  # Design FIR low-pass filter:
  cutoff <- 1 / max(upSamp, downSamp)
  h <- fir1(L - 1, cutoff, window = kaiser(L, beta))
  
  # Function to perform upsampling, filtering, delay compensation, and downsampling.
  upfirdn <- function(xvec, upSamp, downSamp, h) {
    Lx <- length(xvec)
    # Upsample: insert (upSamp-1) zeros between samples.
    x_up <- rep(0, Lx * upSamp)
    x_up[seq(1, Lx * upSamp, by = upSamp)] <- xvec
    
    # Filter the upsampled signal.
    y_filt <- filter(h, 1, x_up)
    
    # Compensate for filter delay: group delay = (length(h)-1)/2.
    delay <- floor((length(h) - 1) / 2)
    if (length(y_filt) > delay) {
      y_filt <- y_filt[(delay + 1):length(y_filt)]
    } else {
      warning("Signal too short for delay compensation. Returning filtered signal without delay compensation.")
    }
    
    # Downsample: take every downSamp-th sample.
    y_down <- y_filt[seq(1, length(y_filt), by = downSamp)]
    return(y_down)
  }
  
  # Apply the upfirdn process to each column of inputMatrix.
  resampled <- apply(inputMatrix, 2, function(col) upfirdn(col, upSamp, downSamp, h))
  
  # If the original input was a vector, convert the output back to a vector.
  if (is.null(original_dim)) {
    resampled <- as.vector(resampled)
  }
  
  return(resampled)
}

# --- getRData Function ---
getRData <- function(data) {
  # data = 1 for original experiment, any other input is fatal  
  if (data == 1) {
    
    # libraries
    library(R.matlab)   # To read .mat data files
    library(mgcv)       # For GAMM analysis
    library(knitr)      # purl(); source .Rmd code for calling functions   
    
    # Base data directory
    base_dir <- "M:/Lab_Shared/Juraj/SpeechEEG_SentencesWithPupillometry/data"
    
    # Miscellaneous directory for other experiment stuff (i.e. eyeUsed table)
    miscDir <- "M:/Lab_Shared/Brady_C/Projects/pupilEEG_fatigue/statistics/scripts/misc"
    
    # Directory to save downsampled data 
    downSampledDataDirectory <- "M:/Lab_Shared/Brady_C/Projects/pupilEEG_fatigue/statistics/scripts/downSampDat"
    
    # Import subject ID function
    purl("M:/Lab_Shared/Brady_C/Projects/pupilEEG_fatigue/statistics/scripts/codebase/getSubjID_RV.Rmd", 
         output = "getSubjID_RV.R")
    source("getSubjID_RV.R")
    
    # Get subject IDs 
    subjIDS <- getSubjList(1)
    
    # Allocate space for data frame for efficiency & a little speed... 
    subject_data <- vector("list", length(subjIDS))
    names(subject_data) <- subjIDS
    
    # Prepare bookkeeping table for eye data selected in saving 
    eyeUsed_Data <- data.frame(
      subject = subjIDS,
      eyeUsed = rep(NA, length(subjIDS)),
      stringsAsFactors = FALSE
    )
    
    # Begin subject loop 
    for (subj in subjIDS) {
      # Construct file path for this subject.
      file_path <- file.path(base_dir, subj, "eyeTrack", "processed", paste0(subj, "_preprocPupilData.mat"))
      
      # Existence check 
      if (!file.exists(file_path)) {
        warning("File does not exist for subject: ", subj)
        stop()  # Fatal error if data does not exist
      }
      
      # Read .mat data file (seems to be time heavy...)
      mat_data <- readMat(file_path)
      
      allDat <- mat_data$allDat
      message("Original allDat dimensions: ", paste(dim(allDat), collapse = " x "))
      current_points <- dim(allDat)[3]
      # Logical check: maximum trial length is 20 seconds; 300Hz sampling should yield 6000 samples.
      if (current_points != 6000) {
        ratio <- 6000 / current_points
        if (ratio %% 1 != 0) {
          warning("Upsampling factor not an int (", ratio, "). Rounding to the nearest int")
          ratio <- round(ratio)
        }
        message("Upsampling allDat from ", current_points, " to 6000 timepoints using upSamp = ", ratio)
        # Create an empty array to hold the upsampled data.
        new_allDat <- array(NA, dim = c(dim(allDat)[1], dim(allDat)[2], 6000))
        # Loop over each trial (first dimension) and each eye channel (second dimension)
        for (i in 1:dim(allDat)[1]) {
          for (j in 1:dim(allDat)[2]) {
            new_allDat[i, j, ] <- resample_signal(as.numeric(allDat[i, j, ]), upSamp = ratio, downSamp = 1, filt = 10, beta = 5)
          }
        }
        allDat <- new_allDat
        message("New allDat dimensions: ", paste(dim(allDat), collapse = " x "))
      }
      
      # Check cond order exists 
      if ("condOrder" %in% names(mat_data)) {
        condOrder <- mat_data$condOrder
      } else {
        warning(subj, ": 'condOrder' not found in .mat file.")
        stop()  # Fatal error; stop
      }
      
      # Check allDat exists 
      if (!"allDat" %in% names(mat_data)) {
        warning(subj, ": 'allDat' not found in the .mat file.")
        stop()
      }
      
      # Compute mean validity for selecting the better eye.
      validityVals <- colMeans(mat_data$allValidityPercentage, na.rm = TRUE)
      
      # Select the better eye based on validity values.
      if (validityVals[1] > validityVals[2]){
        pupilData <- allDat[,1,] 
        eyeUsed <- "L"
      } else {
        pupilData <- allDat[,2,] 
        eyeUsed <- "R"
      }    
      eyeUsed_Data[eyeUsed_Data$subject == subj, "eyeUsed"] <- eyeUsed
      
      # --- New Resampling Step ---
      # Downsample pupil data from 300Hz to 150Hz using the resample_signal function.
      # pupilData is a matrix of [trials x timepoints] (200 x 6000) and we want to resample along time.
      pupilData_down <- t(resample_signal(t(pupilData), upSamp = 1, downSamp = 2, filt = 10, beta = 5))
      
      # Verify dimensions: Should be 200 x 3000 after downsampling.
      if (!all(dim(pupilData_down) == c(200, 3000))) {
        warning("Subject ", subj, ": downsampled data dims != [200 x 3000]")
      }
      
      # Replace original pupilData with the downsampled version.
      pupilData <- pupilData_down
      
      # Check condOrder validity.
      if (is.null(mat_data$condOrder) || any(is.na(mat_data$condOrder))) {
        warning("Subject ", subj, ": condOrder not found or contains NA/NaN values.")
        stop()
      } else {
        condOrder <- mat_data$condOrder
      }
      
      # Save subject's data (pupilData and condOrder) in the list.
      subject_data[[subj]] <- list(pupilData = pupilData, condOrder = condOrder)
      
      # Create subject folder in the downsampled data directory, if it does not exist.
      subjectFolder <- file.path(downSampledDataDirectory, subj)
      if (!dir.exists(subjectFolder)) {
        dir.create(subjectFolder, recursive = TRUE)
      }
      
      today_str <- format(Sys.Date(), "%Y%m%d")
      fileName <- file.path(subjectFolder, paste0("downsampled_", today_str, "_150Hz.rds"))
      saveRDS(pupilData, file = fileName)
      message("Saved downsampled data for ", subj, " at ", fileName)
    }
    
    # Save the eyeUsed_Data table to the misc directory.
    eyeUsedFile <- file.path(miscDir, paste0("eyeUsed_Data_", format(Sys.Date(), "%Y%m%d"), ".csv"))
    write.csv(eyeUsed_Data, eyeUsedFile, row.names = FALSE)
    message("Saved eyeUsed data at ", eyeUsedFile)
    
    # Print dimensions of downsampled subject data.
    for (subj in names(subject_data)) {
      dims <- paste(dim(subject_data[[subj]]$pupilData), collapse = " x ")
      message("Subject: ", subj, " - Downsampled data dimensions: ", dims)
    }
    
    # Return list with pupil data and condition order.
    return(subject_data)
  } else {
    stop("Data value not recognized. Please provide a valid data input.")
  }
}
