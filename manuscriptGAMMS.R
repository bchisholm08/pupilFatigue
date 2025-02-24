# get some libraries that can help with MATLAB 
library(R.matlab)   # To read .mat files
library(mgcv)       # For GAMM analysis

# Set the directory containing your .mat files
mat_dir <- "path/to/your/mat/files"  # change to your actual directory
mat_files <- list.files(mat_dir, pattern = "\\.mat$", full.names = TRUE)

# store model results in a list 
gamm_results <- list()

# Loop over each subject's .mat file
for(file in mat_files) {
  # Read the .mat file
  mat_data <- readMat(file) 
  
  # Extract your data from the .mat file.
  # This example assumes that your .mat file contains a variable named "data"
  # which is organized with columns such as time, pupil size, condition, and trial.
  # Adjust the extraction below to match your file's structure.
  subject_df <- as.data.frame(mat_data$data)
  
  # (Optional) Rename columns if necessary
  # For example, if your columns are in order: time, pupil, condition, trial:
  # names(subject_df) <- c("time", "pupil", "condition", "trial")
  
  # Fit a GAMM model:
  # Here, we model pupil size as a function of a smooth term over time,
  # a parametric effect of condition, and a smooth-by-condition interaction.
  # A random intercept for "trial" is added.
  model <- gamm(pupil ~ s(time) + condition + s(time, by = condition),
                random = list(trial = ~1),
                data = subject_df)
  
  # Use the file name (without extension) as the subject ID
  subject_id <- tools::file_path_sans_ext(basename(file))
  gamm_results[[subject_id]] <- model
}

# To review the summary of each subject's model, you can use:
lapply(gamm_results, summary)
