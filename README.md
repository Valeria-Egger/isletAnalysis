# IsletAnalysisPipeline
## Overview
This project will host the currently under development pipeline for an experimental analysis of Calcium imaging of endocrine cells in the islet (in Danio Rerio). Right now this project includes logic to compute peak detection, extraction of features of detected peaks as well as peak counting on traces. Furthermore some experimental logic on oscillation detection (by computation of Autocorrelation and Morlet Wavelet transformation) is implemented. For usability a script for min-max scaling is provided to run before the detection pipeline. \
**On AI usage** \
AI was used for synthax and debugging help in order to learn the language faster because Matlab is still a new language to me. 

## Pipeline Steps
1. Run "Calculations_CalciumImaging.m" on your data
   1. select the csv file that was provided by ImageJ after Multimeasure of your ROIs
   2. in the Folder that this file is located a folder is going to be created called "scaled"
   3. "scaled" will later have your scaled data
2. Go to mainFunction.m and run this file
3. You may want to adjust the following parameters for the functions
   1. lag - determines how many datapoints are going to be taken into account for the rolling baseline
   2. threshold - determines how many Z-scores over your baseline a signal is counted as signal
   3. influence - how much influence a new datapoint has to the baseline
   4. Minimum signal - a simple threshold which determines a signal under which nothing is going to be counted as a peak no matter if the threshold would allow it or not
   5. exampleCell - let's you plot the peak detection, autocorrelation and wavelet for a cell that you decide for visual inspection

## Input and Output formats
All inputs should be csv files, if you want to run your own scaling over your data before detecting peaks please be sure to arrange the csv in a way that match the output csv of "Calculations_CalciumImaging.m". 

## Dependencies
Matlab R2025b is needed to run this code. Other versions have not been tested.

## Limitations
For now baseline drift is not handled, this will be adressed in the future. However keep in mind that noisy cells or cells with high baseline drift may produce inconsistent results. To adress this visual inspection, parameter tweaking or the Minimum_signal parameter might be helpful. Otherwise consider excluding the cell of the analysis. 
### wavelet function
The frequency limits of the Morlet wavelet transform is determined by the number of timepoints in the dataset and the sampling rate. For the example dataset this would be N (number of timepoints) = 1080 and fs (sampling rate) = 1 (for 1 second between images). 
This leads to frequency limits of 2.3 seconds and 5.5 minutes as min and max respectively. 
To find out your specific frequency limits for your dataset run following command in the command line:
"[minfreq, maxfreq] = cwtfreqbounds(N, fs)"
To convert it to oscillation period simple divide 1 with your limits.
For example in case of the example data (N = 1080 and fs = 1) the command gives minfreq = 0.0031, the maxfreq = 0.4341. This is then converted to: 1/0.0031 = 322.6 seconds and 1/0.4341 = 2.3 seconds. Therefore the wavelet transform for this dataset can reasonable detect oscillations with a period between 2.3 seconds and 5.4 minutes. This was determined reasonable for islet physiology. Please confirm beforehand that your sampling parameters will produce reasonable limits for your application. 
### count Peaks per window
This function will undergo change to allow the user to choose their own time windows rather than relying on the hardcoded provided for this. Right now the hardcoded version supports 1080 timesteps with change (perfusion) every 180 seconds (timepoints). If this is suitable for you at this moment you can simply outcomment the sections to calculate peaks per window and plot it in Matlab. You are also welcome to change this version as needed. 


## Contact
This is just the minimal version of the analysis pipeline. So far it runs and provides the expected output however output is for now only variables and plots and will later be adjusted to csv files that can be exported into R. 
If there are any questions please feel free to reach out.


## Future plans
1. following functions may be included in later versions:
    1. AUC calculations, time to peak
    2. cross-correlation/synchronicity computations
    3. ...
2. in future versions the plotting will be relayed to R with it's own pipeline working with exported csv files of this MatLab pipeline. This ensures clear separation of mathematical computation and visualization and uses the strengths of both programming languages.
3. for future versions this pipeline might be rewritten in python to avoid licenses
4. the provided functions here are stable but may still be subject to change

This project is under active development. Feedback, suggestions and contributions are welcome.

