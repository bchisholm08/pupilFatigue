#-------------------------------------------------------
# Author: Brady M. Chisholm
# University of Minnesota Twin Cities, Dpt. of Psychology
# 3.17.2025
#-------------------------------------------------------

# Load required package
if (!require(signal)) {
  install.packages("signal")
  library(signal)
}

# Debugged version of resample_signal function
resample_signal <- function(inputMatrix, upSamp, downSamp, filt = 10, beta = 5) { 
  # Fix: Added default value for 'filt' to match documentation (default 10)
  # x: Numeric vector or matrix (each column treated as a channel) <inputMatrix>
  # upSamp: Upsampling factor <upSamp>
  # downSamp: Downsampling factor <downSamp>
  # filt: Filter design parameter (default 10) <filt>
  # beta: Kaiser window beta parameter (default 5) <beta = 5>
  
  # If inputMatrix is a vector, convert to a one-column matrix for uniform handling.
  original_dim <- dim(inputMatrix)
  if (is.null(original_dim)) {
    inputMatrix <- matrix(inputMatrix, ncol = 1)  # Fix: Reassign to inputMatrix to ensure uniform handling
  }
  
  # Determine filter length:
  factor <- max(1, downSamp / upSamp)
  L <- as.integer(2 * filt * factor + 1)
  
  # Design FIR low-pass filter:
  cutoff <- 1 / max(upSamp, downSamp)  # Normalized cutoff frequency ensures proper scaling for the resampling
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
      # Fix: Added warning to handle cases where the filtered signal is too short for delay removal.
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

