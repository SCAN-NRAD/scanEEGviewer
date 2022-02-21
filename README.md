# scanEEGviewer
A MATLAB based reader for EEG data which allows simple annotation of epileptiform events for further processing.
- developed on MATLAB 2019b, tested on 2015b, does not run on <= 2014a
- optimized for screen resolution 1920x1200

![alt text](https://raw.githubusercontent.com/SCAN-NRAD/scanEEGviewer/main/ExampPrtScr.png)


## Quick start:
- download script *scanEEGviewer.m*
- download and extract *exampleData.zip* to subfolders *./Data* and *./Data/annot* respectively
- in MATLAB run script *scanEEGviewer.m*


### Usage with EDF:
- instructions can be found in the script *convertEDFtoMAT.m*


## Full data set:
All EEG and annotation data used in our publication[^1] *More than spikes: on the added value of nonlinear iEEG analysis for surgery planning in temporal lobe epilepsy* can be found at http://mia-software.artorg.unibe.ch/scanEEGdata/.

[^1]: Müller M, Dekkers M, Wiest R, Schindler K, Rummel C. More Than Spikes: On the Added Value of Non-linear Intracranial EEG Analysis for Surgery Planning in Temporal Lobe Epilepsy. Front Neurol. 2022 Jan 13;12:741450. doi: 10.3389/fneur.2021.741450.
