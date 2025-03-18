#-------------------------------------------------------
# Author: Brady M. Chisholm
# University of Minnesota Twin Cities, Dpt. of Psychology
# Date: 3.17.2025
#
# Description: This script loads raw MATLAB pupil data, ensures that
#              the data are uniformly sampled at 300Hz (6000 samples for
#              20-second trials), downsamples to 150Hz (3000 samples), 
#              selects the better eye, and saves the processed data.
#-------------------------------------------------------

getRData <- function(data) {
  # data = 1 for original experiment; any other input is fatal  
  if (data == 1) {
    
    # libraries
    library(R.matlab)   # To read .mat data files
    library(mgcv)       # For GAMM analysis
    library(knitr)      # purl(); source .Rmd code for calling functions   
    
    # Base data directory
    base_dir <- "M:/Lab_Shared/Juraj/SpeechEEG_SentencesWithPupillometry/data"
    
    # Directory for miscellaneous output (e.g., eyeUsed table)
    miscDir <- "M:/Lab_Shared/Brady_C/Projects/pupilEEG_fatigue/statistics/scripts/misc"
    
    # Directory to save downsampled data 
    downSampledDataDirectory <- "M:/Lab_Shared/Brady_C/Projects/pupilEEG_fatigue/statistics/scripts/downSampDat"
    
    # Import subject ID function (assumes getSubjID_RV.R has been properly documented)
    purl("M:/Lab_Shared/Brady_C/Projects/pupilEEG_fatigue/statistics/scripts/codebase/getSubjID_RV.Rmd", 
         output = "getSubjID_RV.R")
    source("getSubjID_RV.R")
    
    # Get subject IDs 
    subjIDS <- getSubjList(1)
    
    # Allocate space for data (as a list, with subject IDs as names)
    subject_data <- vector("list", length(subjIDS))
    names(subject_data) <- subjIDS
    
    # Prepare bookkeeping table for the eye selected
    eyeUsed_Data <- data.frame(
      subject = subjIDS,
      eyeUsed = rep(NA, length(subjIDS)),
      stringsAsFactors = FALSE
    )
    
    # --- Begin subject loop ---
    for (subj in subjIDS) {
      # Construct the file path for this subject’s MATLAB file.
      file_path <- file.path(base_dir, subj, "eyeTrack", "processed", paste0(subj, "_preprocPupilData.mat"))
      
      # Existence check (fatal error if file not found)
      if (!file.exists(file_path)) {
        warning("File does not exist for subject: ", subj)
        stop()  # Terminate execution if file is missing
      }
      
      # Read the .mat file (this can be time heavy)
      mat_data <- readMat(file_path)
      
      # Extract the raw data array
      allDat <- mat_data$allDat
      message("Original allDat dimensions: ", paste(dim(allDat), collapse = " x "))
      
      # Check if data are at 300Hz (i.e., 6000 timepoints for 20 sec trials)
      current_points <- dim(allDat)[3]
      if (current_points != 6000) {
        # Calculate the upsampling ratio needed to reach 6000 samples
        ratio <- 6000 / current_points
        if (ratio %% 1 != 0) {
          warning("Upsampling factor is not an integer (", ratio, "). Rounding to the nearest integer.")
          ratio <- round(ratio)
        }
        message("Upsampling allDat from ", current_points, " to 6000 timepoints using upSamp = ", ratio)
        
        # Create an empty array for the upsampled data
        new_allDat <- array(NA, dim = c(dim(allDat)[1], dim(allDat)[2], 6000))
        # Loop over trials (first dim) and channels (second dim) to apply inline resampling
        for (i in 1:dim(allDat)[1]) {
          for (j in 1:dim(allDat)[2]) {
            # --- Inline Resampling Logic ---
            # Extract the original time series for this trial/channel:
            xvec <- as.numeric(allDat[i, j, ])
            Lx <- length(xvec)
            # Upsample: create a vector with zeros inserted between samples.
            x_up <- rep(0, Lx * ratio)
            x_up[seq(1, Lx * ratio, by = ratio)] <- xvec
            
            # Determine filter length.
            # For upSamp = ratio and downSamp = 1, factor = max(1, 1/ratio) will be 1 if ratio > 1.
            L <- as.integer(2 * 10 * 1 + 1)  # Using filt = 10
            # Design an FIR low-pass filter with a Kaiser window.
            cutoff <- 1 / max(ratio, 1)       # Normalized cutoff frequency
            h <- fir1(L - 1, cutoff, window = kaiser(L, 5))  # beta = 5
            
            # Filter the upsampled signal.
            y_filt <- filter(h, 1, x_up)
            
            # Compensate for filter delay.
            delay <- floor((length(h) - 1) / 2)
            if (length(y_filt) > delay) {
              y_filt <- y_filt[(delay + 1):length(y_filt)]
            } else {
              warning("Signal too short for delay compensation. Returning filtered signal without delay removal.")
            }
            
            # With downSamp = 1, the result is y_filt.
            # Ensure exactly 6000 samples (truncate if necessary)
            if (length(y_filt) >= 6000) {
              y_down <- y_filt[1:6000]
            } else {
              # If fewer than 6000 samples, pad with NA values.
              y_down <- c(y_filt, rep(NA, 6000 - length(y_filt)))
            }
            # Save the processed time series.
            new_allDat[i, j, ] <- y_down
          }
        }
        allDat <- new_allDat
        message("New allDat dimensions: ", paste(dim(allDat), collapse = " x "))
      }
      
      # Ensure the 'condOrder' field exists
      if ("condOrder" %in% names(mat_data)) {
        condOrder <- mat_data$condOrder
      } else {
        warning(subj, ": 'condOrder' not found in .mat file.")
        stop()
      }
      
      # Compute mean validity to choose the best eye.
      validityVals <- colMeans(mat_data$allValidityPercentage, na.rm = TRUE)
      if (validityVals[1] > validityVals[2]) {
        pupilData <- allDat[, 1, ]
        eyeUsed <- "L"
      } else {
        pupilData <- allDat[, 2, ]
        eyeUsed <- "R"
      }
      eyeUsed_Data[eyeUsed_Data$subject == subj, "eyeUsed"] <- eyeUsed
      
      # --- Downsampling Step ---
      # Now that pupilData is a matrix of [trials x timepoints] (200 x 6000),
      # downsample from 300Hz to 150Hz by taking every 2nd sample.
      # We implement the same filtering/delay compensation approach inline via transposition.
      # For clarity, here we use a simple approach with our inline resampling logic:
      # Transpose pupilData so that time points are in columns:
      pupilData_T <- t(pupilData)
      # Downsample by selecting every 2nd sample:
      pupilData_down <- pupilData_T[seq(1, nrow(pupilData_T), by = 2), ]
      # Transpose back so the dimensions are [trials x timepoints]:
      pupilData_down <- t(pupilData_down)
      
      if (!all(dim(pupilData_down) == c(200, 3000))) {
        warning("Subject ", subj, ": downsampled data dims != [200 x 3000]")
      }
      pupilData <- pupilData_down
      
      # Final check on condOrder
      if (is.null(mat_data$condOrder) || any(is.na(mat_data$condOrder))) {
        warning("Subject ", subj, ": condOrder not found or contains NA/NaN values.")
        stop()
      } else {
        condOrder <- mat_data$condOrder
      }
      
      # Save the subject's processed data.
      subject_data[[subj]] <- list(pupilData = pupilData, condOrder = condOrder)
      
      # Create folder for subject if it doesn't exist.
      subjectFolder <- file.path(downSampledDataDirectory, subj)
      if (!dir.exists(subjectFolder)) {
        dir.create(subjectFolder, recursive = TRUE)
      }
      
      today_str <- format(Sys.Date(), "%Y%m%d")
      fileName <- file.path(subjectFolder, paste0("downsampled_", today_str, "_150Hz.rds"))
      saveRDS(pupilData, file = fileName)
      message("Saved downsampled data for ", subj, " at ", fileName)
    }
    
    # Save the eyeUsed bookkeeping table.
    eyeUsedFile <- file.path(miscDir, paste0("eyeUsed_Data_", format(Sys.Date(), "%Y%m%d"), ".csv"))
    write.csv(eyeUsed_Data, eyeUsedFile, row.names = FALSE)
    message("Saved eyeUsed data at ", eyeUsedFile)
    
    # Print dimensions for each subject's downsampled data.
    for (subj in names(subject_data)) {
      dims <- paste(dim(subject_data[[subj]]$pupilData), collapse = " x ")
      message("Subject: ", subj, " - Downsampled data dimensions: ", dims)
    }
    
    # Return the list with pupil data and condition order.
    return(subject_data)
    
  } else {
    stop("Data value not recognized. Please provide a valid data input.")
  }
}
