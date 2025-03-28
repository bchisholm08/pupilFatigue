# pupilFatigue Analysis

This repository is dedicated to statistics from the pupillary fatigue project in the APC Lab at UMN-TC. 

Utilizes `mgcv` R package for facilitating GAMM's, and `R.matlab` for handling data files. 

Reach out to [chish071@umn.edu](mailto:chish071@umn.edu) with any questions!

### Code Guide

>`gammManuscriptStats.Rmd`

Loads needed libraries, then uses two external scripts, `getRData()` and `resample_signal()`.
After this, subject data is extracted for statistics to be conducted. 

>`resampleSignal.Rmd`

Changes sample rate of given data. Accepts either a numeric vector or a matrix as input. If given a vector, data is converted 
to a one-column matrix. 

The filter (`L`) is based on the `filt` parameter of the function and the ratio between the downsampling and upsampling function inputs. 

A finite impulse response (FIR) low-pass filter is designed wiith `fir1` using a Kaiser window. 
Cutoff frequency is determined by the `max` of the up/downsample inputs. This ratio tells us what our resulting data will be sampled at. 

Within this function, a helper function called `upfirdn`:

- ***Upsamples*** by inserting zeros between original data samples 
- ***Filters*** with the designated filter on the upsamples signal 
- ***Computes Delay Compensation*** by computing half of the filter length minus one, removing delay introduced into the signal.
- ***Downsamples*** by selecting every `downSamp`-th sample from filtered signal

This process is applied to each column of the input matrix, which here represents *INSERT*. If the original data was a vector, a vector is output.

>`getSubjID_RV.Rmd`

Gathers the list of subject ID's associated with experiment participants. 

Vector of subject ID's must be manually updated. Used to iterate subjects in the data process.

>`getData_RV.Rmd`

Handles loading, processing, and saving of pupillometry data with a subject-wise loop.

Gathers subject ID list, data file pathways, and pathways to save to. 

After set up, the subject-wise loop builds a `MATLAB` file pathway for each 
participant and reads their `.mat` data file, which is read into R with `readMat()`

After `.mat` data is read, the quality of the data is checked and data from one eye 
is selected for use in statistics, based on a better average validity percentage over 
over the course of the experiment. 

After better data is gathered for the subject, data is transposed and resampled. 

`resample_signal` is used with an upsampling input of 1, and downsampling input of 2, which reduces sampling rate by 1/2. 
This changes sampling from 300Hz to 150Hz. A final check confirms matrix dimensions against expected dimensions.

Data is saved int the `subject_data` list, along with the `condOrder` variable. Data is saved into subject-specific directories and downsampled data 
is stored in `.rds` file format.  

A table, `eyeUsedData`, is saved for tracking which eye data was used for each subject. This is a two column table of subject ID, and L/R, indicating which eye was used. 