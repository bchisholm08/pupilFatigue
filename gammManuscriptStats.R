#-------------------------------------------------------
# Author: Brady M. Chisholm
# University of Minnesota Twin Cities, Dpt. of Psychology
# Date: 3.17.2025
#
# Description: This script loads processed pupil data using the 
#              integrated getRData() function (which now includes the
#              resample_signal function), combines the data into a 
#              long-format data frame, computes summary statistics, 
#              and produces several exploratory plots.
#-------------------------------------------------------

library(mgcv)       # For GAMM analysis 
library(dplyr)      # Data manipulations & formatting
library(tidyr)      # Data manipulations & formatting
library(ggplot2)    # Figures

# Source the getData_RV.R script, which now integrates the resample_signal function
purl("M:/Lab_Shared/Brady_C/Projects/pupilEEG_fatigue/statistics/scripts/codebase/getSubjID_RV.R", 
     output = "getSubjID_RV.R")
source("getSubjID_RV.R")
# No need to source "resampleSignal.R" separately

# Load processed subject data
subjDataList <- getSubjID_RV(1)
# Data dimension assertions are handled within getRData()
summary(subjDataList)

# Combine the downsampled data from all subjects into one long-format data frame
all_data <- lapply(names(subjDataList), function(subj) {
  data_subj <- subjDataList[[subj]]
  
  # Extract the pupil data (trials x timepoints) and condition order (vector)
  pupilData <- data_subj$pupilData    
  condOrder <- data_subj$condOrder      
  
  # Convert the pupil matrix to a data frame; using cbind() avoids issues with vector recycling
  df <- as.data.frame(pupilData)
  df <- cbind(df, trial = 1:nrow(df), subject = subj, condOrder = condOrder)
  
  # Reshape from wide (each column = a time point, e.g., V1, V2, ...) to long format
  df_long <- pivot_longer(df,
                          cols = starts_with("V"),
                          names_to = "time",
                          values_to = "pupil")
  
  # Convert the time variable (e.g., "V1") to a numeric index
  df_long$time <- as.numeric(gsub("V", "", df_long$time))
  
  return(df_long)
})

# Bind all subject data together
combined_data <- bind_rows(all_data)

# Inspect structure and summary of the combined data
str(combined_data)
summary(combined_data)

# Save the combined data for future use
saveRDS(combined_data, file = "combined_pupil_data.rds")

# --- Summary Statistics ---

# Compute summary stats for each subject/trial (mean and standard deviation)
trial_stats <- combined_data %>%
  group_by(subject, trial) %>%
  summarize(mean_pupil = mean(pupil, na.rm = TRUE),
            sd_pupil   = sd(pupil, na.rm = TRUE),
            .groups = "drop")

# Compute overall summary stats across all subjects and trials
overall_stats <- combined_data %>%
  summarize(mean_pupil = mean(pupil, na.rm = TRUE),
            sd_pupil   = sd(pupil, na.rm = TRUE),
            min_pupil  = min(pupil, na.rm = TRUE),
            max_pupil  = max(pupil, na.rm = TRUE))

print(trial_stats)
print(overall_stats)

# --- Plots ---

# Mean pupil time series per subject
subject_time_series <- combined_data %>%
  group_by(subject, time) %>%
  summarize(mean_pupil = mean(pupil, na.rm = TRUE), .groups = "drop")

ggplot(subject_time_series, aes(x = time, y = mean_pupil, group = subject)) +
  geom_line(alpha = 0.5) +
  labs(title = "Mean Pupil Time Series per Subject",
       x = "Time (sample index)",
       y = "Mean Pupil Diameter")

# Distribution of pupil measurements
ggplot(combined_data, aes(x = pupil)) +
  geom_histogram(bins = 50, fill = "grey", color = "black") +
  labs(title = "Distribution of Pupil Measurements",
       x = "Pupil Diameter",
       y = "Frequency")

# Mean pupil time series by condition
condition_time_series <- combined_data %>%
  group_by(condOrder, time) %>%
  summarize(mean_pupil = mean(pupil, na.rm = TRUE), .groups = "drop")

ggplot(condition_time_series, aes(x = time, y = mean_pupil, color = as.factor(condOrder))) +
  geom_line() +
  labs(title = "Mean Pupil Time Series by Condition",
       x = "Time (sample index)",
       y = "Mean Pupil Diameter",
       color = "Condition")
