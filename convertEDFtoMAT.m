
%%%
% this script reads an EDF file and generates the variables necessary
% to use the data in our EEG-viewer scanEEGviewer.m
% (https://github.com/SCAN-NRAD/scanEEGviewer)
% 
% to read the EDF data, a script by Andrey Shapkin can be used, available at:
% https://www.mathworks.com/matlabcentral/fileexchange/38641-reading-and-saving-of-data-in-the-edf
% Andrey Shapkin (2021). Reading and saving of data in the EDF+
% MATLAB Central File Exchange. Retrieved July 28, 2021. 
% 
% - only channels with the same sampling frequency as the first channel are included
% - further channels containing other data than EEG (e.g. ECG) may also be excluded before use
%
% (c) Michael Müller,	July 2021
%	ORCiD: 0000-0002-6915-4820
%%%


%%% read data
[EEG, header] = readEDF('./some-EDF-file.edf');

%%% select and convert data
fs = header.samplerate(1);
selChIdx = header.samplerate == fs;
EEG = [EEG{selChIdx}];
ELabel = header.labels(selChIdx);

%%% use data
scanEEGviewer(EEG, fs, ELabel)


