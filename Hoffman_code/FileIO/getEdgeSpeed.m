function [edgeSpeed] = getEdgeSpeed()
% This function is used to get edge speed stored by
% leadingEdgeAnalysis.m

% Locate files
myDirectory = uigetdir();
myFile = dir(fullfile(myDirectory,'*.mat'));
myFile = myFile(1).name;
fitLength = 0;
load(fullfile(myDirectory,myFile),'dA_orig','fitLength');

% Convert to um/hr
edgeSpeed = dA_orig./fitLength(1:end-1);

end