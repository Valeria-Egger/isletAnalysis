# IsletAnalysisPipeline
## Overview
This project will host the currently under development pipeline for an experimental analysis of Calcium imaging of endocrine cells in the islet (in Danio Rerio). Right now this project includes logic to compute peak detection, extraction of features of detected peaks as well as peak counting on traces. Furthermore some experimental logic on oscillation detection (by computation of Autocorrelation and Morlet Wavelet transformation) is implemented. For usability a script for min-max scaling is provided to run before the detection pipeline.

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

