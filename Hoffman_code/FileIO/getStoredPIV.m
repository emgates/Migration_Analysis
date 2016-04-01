function [PIV_vectors, PIV_bins] = getStoredPIV(binRange)
% This function is used to generate PIV vectors stored by
% leadingEdgeAnalysis.m

% Locate files
myDirectory = uigetdir();
myFile = dir(fullfile(myDirectory,'*.mat'));
myFile = myFile(1).name;
datetime = 0;
load(fullfile(myDirectory,myFile),'PIV_vectors','datetime');
t = datetime(2:end)-datetime(1:end-1);

% Convert to um/frame
PIV_vectors(:,1:2) = PIV_vectors(:,1:2)*0.32;

% Convert to um/hr
for k = 1:size(PIV_vectors,1)
    PIV_vectors(k,1:2) = PIV_vectors(k,1:2)/t(PIV_vectors(k,5));
end

PIV_bins = histc(sqrt(PIV_vectors(:,1).^2+PIV_vectors(:,2).^2),binRange);
PIV_bins = PIV_bins/sum(PIV_bins);

end